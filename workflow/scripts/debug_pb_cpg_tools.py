

import altair as alt
import polars as pl
import numpy as np
import pandas as pd

pl.Config.set_tbl_rows(20)
pl.Config.set_tbl_cols(20)

COLS = [
    "chromosome",
    "position",
    "varlo_0.01_methylation_rep1",
    "varlo_0.01_methylation_rep2",
    "sample",
]


def load(path, source, extra_cols=None):
    """Read a parquet file, tag it with its source, rename varlo columns."""
    df = pl.read_parquet(path)
    # Do not use loci with coverage < 30
    df = df.with_columns(
        pl.when(
            pl.col("varlo_0.01_format_rep1").str.split(":").list.get(0).cast(pl.Int64)
            < 30
        )
        .then(pl.lit(None))
        .otherwise(pl.col("varlo_0.01_methylation_rep1"))
        .alias("varlo_0.01_methylation_rep1"),
        pl.when(
            pl.col("varlo_0.01_format_rep2").str.split(":").list.get(0).cast(pl.Int64)
            < 30
        )
        .then(pl.lit(None))
        .otherwise(pl.col("varlo_0.01_methylation_rep2"))
        .alias("varlo_0.01_methylation_rep2"),
    )
    df = df.select(*COLS, *(extra_cols or []))
    if "sample" not in df.columns or source != "illumina":
        df = df.with_columns(pl.lit(source).alias("sample"))
    return df.rename(
        {
            "varlo_0.01_methylation_rep1": f"varlo_{source}_rep1",
            "varlo_0.01_methylation_rep2": f"varlo_{source}_rep2",
        }
    )


def to_long(df, source):
    """Convert a wide per-sample dataframe to long format (one row per replicate)."""
    return df.select(
        "chromosome",
        "position",
        "sample",
        "meth_type",
        pl.col(f"varlo_{source}_rep1").alias("rep1"),
        pl.col(f"varlo_{source}_rep2").alias("rep2"),
    ).unpivot(
        index=["chromosome", "position", "sample", "meth_type"],
        on=["rep1", "rep2"],
        variable_name="replicate",
        value_name="methylation",
    )


def to_long_pacbio(df):
    varlo_long = (
        df.select(
            "chromosome",
            "position",
            "meth_type",
            pl.col("varlo_pacbio_rep1").alias("rep1"),
            pl.col("varlo_pacbio_rep2").alias("rep2"),
        )
        .with_columns(pl.lit("pacbio-varlo").alias("sample"))
        .unpivot(
            index=["chromosome", "position", "sample", "meth_type"],
            on=["rep1", "rep2"],
            variable_name="replicate",
            value_name="methylation",
        )
    )
    pbcpg_long = (
        df.select(
            "chromosome",
            "position",
            "meth_type",
            pl.col("pb_CpG_tools_methylation_rep1").alias("rep1"),
            pl.col("pb_CpG_tools_methylation_rep2").alias("rep2"),
        )
        .with_columns(pl.lit("pacbio-pb_cpg").alias("sample"))
        .unpivot(
            index=["chromosome", "position", "sample", "meth_type"],
            on=["rep1", "rep2"],
            variable_name="replicate",
            value_name="methylation",
        )
    )
    return pl.concat([varlo_long, pbcpg_long])

def find_interesting_positions(df, condition_pb, condition_varlo):
    df = df.select(
        "chromosome",
        "position",
        pl.col("pb_CpG_tools_methylation_rep1").alias("pb_cpg_rep1"),
        pl.col("pb_CpG_tools_methylation_rep2").alias("pb_cpg_rep2"),
        pl.col("varlo_pacbio_rep1").alias("varlo_rep1"),
        pl.col("varlo_pacbio_rep2").alias("varlo_rep2"),
    )

    interesting_positions = (
        df
        .unpivot(
            index=["chromosome", "position"],
            on=["pb_cpg_rep1", "pb_cpg_rep2", "varlo_rep1", "varlo_rep2"],
            variable_name="variable",
            value_name="value",
        )
        .with_columns(
            pl.col("variable").str.extract(r"(rep\d)$").alias("replicate"),
            pl.col("variable").str.extract(r"^(pb_cpg|varlo)").alias("metric"),
        )
        .pivot(on="metric", index=["chromosome", "position", "replicate"], values="value")
        .filter(condition_pb(pl.col("pb_cpg")) & condition_varlo(pl.col("varlo")))
        .select("chromosome", "position", "replicate")
    )
    return interesting_positions.select("chromosome", "position").unique()

def merge_and_filter_dfs(pacbio, nanopore, illumina, low_meth_positions, high_meth_positions, sample_coverage):
    pacbio_low, nanopore_low, illumina_low = (
        df.join(low_meth_positions, on=["chromosome", "position"], how="inner") .with_columns(pl.lit("low").alias("meth_type"))
        for df in (pacbio, nanopore, illumina)
    )
    pacbio_high, nanopore_high, illumina_high = (
        df.join(high_meth_positions, on=["chromosome", "position"], how="inner").with_columns(pl.lit("high").alias("meth_type"))
        for df in (pacbio, nanopore, illumina)
    )

    # Combine into long format
    combined = (
        pl.concat(
            [
                to_long_pacbio(pacbio_low),
                to_long_pacbio(pacbio_high),
                to_long(nanopore_low, "nanopore"),
                to_long(nanopore_high, "nanopore"),
                to_long(illumina_low, "illumina"),
                to_long(illumina_high, "illumina"),
            ]
        ).drop_nulls(subset=["methylation"])
        # .join(low_meth_replicates, on=["chromosome", "position", "replicate"], how="inner")
    )

    # Keep positions with methylation values in at least x of 9 samples
    sample_coverage = (
        combined.group_by(["chromosome", "position"])
        .agg(pl.col("sample").n_unique().alias("n_samples"))
        .filter(pl.col("n_samples") >= sample_coverage)
        .select("chromosome", "position")
    )
    combined = combined.join(sample_coverage, on=["chromosome", "position"], how="inner")
    return combined


def plot_histo(combined, title):
    binned = (
        combined.with_columns(
            ((pl.col("methylation") / 5).round(0) * 5).alias("meth_bin")
        ).group_by(["meth_type", "meth_bin"])
        .agg(pl.len().alias("count"))
        .with_columns(
            (pl.col("count") / pl.col("count").sum().over(["meth_type"]))
            .alias("percentage")
        )
    )
    histo = (
        alt.Chart(binned.to_pandas())
        .mark_bar()
        .encode(
            x=alt.X("meth_bin:O", title="Binned methylation" if title=="Fully Methylated" else None),
            xOffset=alt.XOffset("meth_type:N"),
            y=alt.Y("percentage:Q", title="Percentage", axis=alt.Axis(format="%")),
            color=alt.Color("meth_type:N",  scale=alt.Scale(
                range=["#D81B60", "#1E88E5"],
            ),title="Varlo prediction"),
        )
        .properties(height=400, width=800, title=title)
    )
    # chart = alt.vconcat(strip, histo).resolve_scale(color="independent")
    return histo

# Plot scatter plot with methylation per replicate
def build_scatter_data(pacbio, low_positions, high_positions):
    """Join pacbio on low/high positions, keep rep1/rep2 as cols,
    split into two tool rows (varlo, pb_cpg) per locus."""
    pacbio_low = pacbio.join(
        low_positions, on=["chromosome", "position"], how="inner"
    ).with_columns(pl.lit("low").alias("meth_type"))
    pacbio_high = pacbio.join(
        high_positions, on=["chromosome", "position"], how="inner"
    ).with_columns(pl.lit("high").alias("meth_type"))
    combined = pl.concat([pacbio_low, pacbio_high])

    varlo = combined.select(
        "chromosome",
        "position",
        "meth_type",
        pl.col("varlo_pacbio_rep1").alias("rep1"),
        pl.col("varlo_pacbio_rep2").alias("rep2"),
    ).with_columns(pl.lit("varlo").alias("tool"))

    pbcpg = combined.select(
        "chromosome",
        "position",
        "meth_type",
        pl.col("pb_CpG_tools_methylation_rep1").alias("rep1"),
        pl.col("pb_CpG_tools_methylation_rep2").alias("rep2"),
    ).with_columns(pl.lit("pb_cpg").alias("tool"))

    return pl.concat([varlo, pbcpg]).drop_nulls(subset=["rep1", "rep2"])


def plot_scatter(data, title, pb_filter):
    chart = (
        alt.Chart(data.to_pandas())
        .mark_point(filled=True, size=30, opacity=0.3)
        .encode(
            x=alt.X("rep1:Q", title="Replicate 1" if title == "Fully Methylated" else None),
            y=alt.Y("rep2:Q", title="Replicate 2"),
            color=alt.Color("tool:N", title="Tool", scale=alt.Scale(range=["#8D9279", "#05AA8F"])),
            shape=alt.Shape("meth_type:N", title="Varlo prediction"),
            tooltip=["chromosome", "position", "tool", "meth_type", "rep1", "rep2"],
        )
        .properties(height=400, width=400, title=title)
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

    chart = chart + vline + hline

    return chart

# Load data
pacbio = load(
    snakemake.input.pacbio,
    "pacbio",
    ["pb_CpG_tools_methylation_rep1", "pb_CpG_tools_methylation_rep2"],
)
nanopore = load(snakemake.input.nanopore, "nanopore")
illumina = load(snakemake.input.illumina, "illumina")

# Find positions/replicates where pb_CpG_tools < 5 but varlo > 10
low_meth_positions = find_interesting_positions(pacbio, lambda v: v > 95, lambda v: v <= 90)
high_meth_positions = find_interesting_positions(pacbio, lambda v: v > 95, lambda v: v == 100)

# Find positions/replicates where pb_CpG_tools > 95 but varlo < 10
low_unmeth_positions = find_interesting_positions(pacbio, lambda v: v < 5, lambda v: v == 0)
high_unmeth_positions = find_interesting_positions(pacbio, lambda v: v < 5, lambda v: v >= 10)


# pacbio = pacbio.drop("pb_CpG_tools_methylation_rep1", "pb_CpG_tools_methylation_rep2")
# Restrict to positions of interest
combined_unmeth = merge_and_filter_dfs(pacbio, nanopore, illumina, low_unmeth_positions, high_unmeth_positions, 0)
combined_meth = merge_and_filter_dfs(pacbio, nanopore, illumina, low_meth_positions, high_meth_positions, 0)


histo_unmeth = plot_histo(combined_unmeth, "Unmethylated")
histo_meth = plot_histo(combined_meth, "Fully Methylated")
histo = alt.vconcat(histo_unmeth, histo_meth).resolve_scale(x='shared', y='shared')




scatter_data_unmeth = build_scatter_data(pacbio, low_unmeth_positions, high_unmeth_positions)
scatter_data_meth = build_scatter_data(pacbio, low_meth_positions, high_meth_positions)

scatter_unmeth = plot_scatter(scatter_data_unmeth, "Unmethylated", 5)
scatter_meth = plot_scatter(scatter_data_meth, "Fully Methylated", 95)
scatter = alt.vconcat(scatter_unmeth, scatter_meth).resolve_scale(color="shared", shape="shared")


chart = (
    alt.hconcat(scatter, histo)
    .resolve_scale(color="independent", shape="independent")
    .configure_axis(labelFontSize=18, titleFontSize=20)
    .configure_legend(labelFontSize=18, titleFontSize=20)
    .configure_title(fontSize=22, subtitleFontSize=18)
    .configure_header(labelFontSize=18, titleFontSize=20)
)
chart.save(snakemake.output[0])
