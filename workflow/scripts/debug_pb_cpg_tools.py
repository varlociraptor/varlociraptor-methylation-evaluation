import sys

import altair as alt
import pandas as pd
import polars as pl

sys.stderr = open(snakemake.log[0], "w")

pl.Config.set_tbl_rows(20)
pl.Config.set_tbl_cols(20)

# Same colors as in the other plots
TOOL_COLORS = {"pb-CpG-tools": "#8D9279", "Varlociraptor": "#05AA8F"}

min_coverage = snakemake.params.min_coverage


def read_coverage(coverage_file, rep_name):
    """Read a mosdepth regions file and return per-position coverage for one replicate."""
    return (
        pl.read_csv(
            coverage_file,
            separator="\t",
            has_header=False,
            new_columns=["chromosome", "pos_start", "pos_end", "coverage"],
        )
        .with_columns(
            pl.col("chromosome").cast(pl.String),
            ((pl.col("pos_end") + pl.col("pos_start")) / 2)
            .cast(pl.Int64)
            .alias("position"),
        )
        .select(["chromosome", "position", "coverage"])
        .rename({"coverage": f"coverage_{rep_name}"})
    )


def load_coverage(samples, paths_rep1, paths_rep2):
    """Read the coverage of both replicates for every sample.
    Both files of a sample contain the same candidate positions.
    """
    sample_dfs = []
    for sample, path_coverage_rep1, path_coverage_rep2 in zip(samples, paths_rep1, paths_rep2):
        rep1 = read_coverage(path_coverage_rep1, "rep1")
        rep2 = read_coverage(path_coverage_rep2, "rep2")
        sample_df = rep1.join(rep2, on=["chromosome", "position"]).with_columns(
            pl.lit(sample).alias("sample")
        )
        sample_dfs.append(sample_df)
    return pl.concat(sample_dfs)


def load_platform(platform, callers):
    """Read the methylation calls of a platform, one column per caller and replicate.

    `callers` maps the column name used here to the caller name in the parquet file.
    Values with a mosdepth coverage below min_coverage in their replicate are removed,
    so the filter is the same for all callers and does not depend on Varlociraptor.
    """
    coverage = load_coverage(
        snakemake.params[f"{platform}_samples"],
        snakemake.input[f"{platform}_coverages_rep1"],
        snakemake.input[f"{platform}_coverages_rep2"],
    )
    platform_df = pl.scan_parquet(snakemake.input[platform]).select(
        "chromosome",
        "position",
        "sample",
        *(
            f"{caller}_methylation_{rep}"
            for caller in callers.values()
            for rep in ["rep1", "rep2"]
        ),
    ).collect()
    df = platform_df.join(coverage, on=["sample", "chromosome", "position"], how="left")
    return df.select(
            "chromosome",
            "position",
            "sample",
            *(
                pl.when(pl.col(f"coverage_{rep}") >= min_coverage)
                .then(pl.col(f"{caller}_methylation_{rep}"))
                .alias(f"{name}_{rep}")
                for name, caller in callers.items()
                for rep in ["rep1", "rep2"]
            ),
        )



def find_loci(pacbio, pb_condition, varlo_condition, meth_type):
    """Loci where both conditions hold in the same replicate, in at least one replicate."""
    return (
        pacbio.filter(
            pl.any_horizontal(
                pb_condition(pl.col(f"pb_cpg_{rep}"))
                & varlo_condition(pl.col(f"varlo_{rep}"))
                for rep in ["rep1", "rep2"]
            )
        )
        .select("chromosome", "position")
        .unique()
        .with_columns(pl.lit(meth_type).alias("meth_type"))
    )


def count_subtitle(df, with_values):
    """Number of loci (and values) per label, e.g. 'low: 727 loci (1,881 values)'."""

    def label_text(meth_type):
        subset = df.filter(pl.col("meth_type") == meth_type)
        text = f"{meth_type}: {subset.n_unique(subset=['chromosome', 'position']):,} loci"
        if with_values:
            text += f" ({subset.height:,} values)"
        return text

    low, high = label_text("extreme"), label_text("unextreme")
    return f"{low} · {high}"


def plot_scatter(df, title, pb_filter):
    # One point per locus and tool; loci need a value in both replicates
    data = pl.concat(
        [
            df.select(
                "chromosome",
                "position",
                "meth_type",
                pl.col("varlo_rep1").alias("rep1"),
                pl.col("varlo_rep2").alias("rep2"),
                pl.lit("Varlociraptor").alias("tool"),
            ),
            df.select(
                "chromosome",
                "position",
                "meth_type",
                pl.col("pb_cpg_rep1").alias("rep1"),
                pl.col("pb_cpg_rep2").alias("rep2"),
                pl.lit("pb-CpG-tools").alias("tool"),
            ),
        ]
    ).drop_nulls(subset=["rep1", "rep2"])
    chart = (
        alt.Chart(data.to_pandas())
        .mark_point(filled=True, size=30, opacity=0.3)
        .encode(
            x=alt.X(
                "rep1:Q", title="Replicate 1" if title == "Fully Methylated" else None
            ),
            y=alt.Y("rep2:Q", title="Replicate 2"),
            color=alt.Color(
                "tool:N",
                title="Tool",
                scale=alt.Scale(
                    domain=list(TOOL_COLORS), range=list(TOOL_COLORS.values())
                ),
            ),
            shape=alt.Shape("meth_type:N", title="Varlo prediction"),
            tooltip=["chromosome", "position", "tool", "meth_type", "rep1", "rep2"],
        )
        .properties(
            height=400,
            width=400,
            title=alt.Title(title, subtitle=count_subtitle(data, with_values=False)),
        )
    )

    threshold_df = pd.DataFrame({"val": [pb_filter]})

    vline = (
        alt.Chart(threshold_df)
        .mark_rule(color="red", size=1, strokeDash=[4, 4])
        .encode(x=alt.X("val:Q"))
    )
    hline = (
        alt.Chart(threshold_df)
        .mark_rule(color="red", size=1, strokeDash=[4, 4])
        .encode(y=alt.Y("val:Q"))
    )

    return chart + vline + hline

def plot_histo(df, title):
    # One value per replicate; replicates without a value are not counted
    values = df.unpivot(
        index=["chromosome", "position", "sample", "meth_type"],
        on=["independent_rep1", "independent_rep2"],
        value_name="methylation",
    ).drop_nulls(subset=["methylation"])
    binned = (
        values.with_columns(
            ((pl.col("methylation") / 5).round(0) * 5).alias("meth_bin")
        )
        .group_by(["meth_type", "meth_bin"])
        .agg(pl.len().alias("count"))
        .with_columns(
            (pl.col("count") / pl.col("count").sum().over(["meth_type"])).alias(
                "percentage"
            )
        )
    )
    return (
        alt.Chart(binned.to_pandas())
        .mark_bar()
        .encode(
            x=alt.X(
                "meth_bin:O",
                title="Binned methylation" if title == "Fully Methylated" else None,
            ),
            xOffset=alt.XOffset("meth_type:N"),
            y=alt.Y("percentage:Q", title="Percentage", axis=alt.Axis(format="%")),
            color=alt.Color(
                "meth_type:N",
                scale=alt.Scale(range=["#D81B60", "#1E88E5"]),
                title="Varlo prediction",
            ),
        )
        .properties(
            height=400,
            width=800,
            # Percentages refer to values, i.e. loci x samples x replicates
            title=alt.Title(title, subtitle=count_subtitle(values, with_values=True)),
        )
    )


pacbio = load_platform("pacbio", {"varlo": "varlo_0.01", "pb_cpg": "pb_CpG_tools"})
# Callers independent of Varlociraptor, used to check the PacBio loci on the other platforms.
# Same column names, so both platforms fit into one table.
independent = pl.concat(
    [
        load_platform("nanopore", {"independent": "modkit"}),
        load_platform("illumina", {"independent": "methylDackel"}),
    ]
)

# pb-CpG-tools < 5; low: Varlociraptor = 0, high: Varlociraptor >= 10
unmeth_loci = pl.concat(
    [
        find_loci(pacbio, lambda v: v < 5, lambda v: v == 0, "extreme"),
        find_loci(pacbio, lambda v: v < 5, lambda v: v >= 10, "unextreme"),
    ]
)
# pb-CpG-tools > 95; low: Varlociraptor <= 90, high: Varlociraptor = 100
meth_loci = pl.concat(
    [
        find_loci(pacbio, lambda v: v > 95, lambda v: v <= 90, "unextreme"),
        find_loci(pacbio, lambda v: v > 95, lambda v: v == 100, "extreme"),
    ]
)


scatter = alt.vconcat(
    plot_scatter(pacbio.join(unmeth_loci, on= ["chromosome", "position"]), "Unmethylated", 5),
    plot_scatter(pacbio.join(meth_loci, on= ["chromosome", "position"]), "Fully Methylated", 95),
).resolve_scale(color="shared", shape="shared")
histo = alt.vconcat(
    plot_histo(independent.join(unmeth_loci, on= ["chromosome", "position"]), "Unmethylated"),
    plot_histo(independent.join(meth_loci, on= ["chromosome", "position"]), "Fully Methylated"),
).resolve_scale(x="shared", y="shared")

chart = (
    alt.hconcat(scatter, histo)
    .resolve_scale(color="independent", shape="independent")
    .configure_axis(labelFontSize=18, titleFontSize=20)
    .configure_legend(labelFontSize=18, titleFontSize=20)
    .configure_title(fontSize=22, subtitleFontSize=18)
    .configure_header(labelFontSize=18, titleFontSize=20)
)
chart.save(snakemake.output[0])
