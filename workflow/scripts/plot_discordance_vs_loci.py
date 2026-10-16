import sys

import altair as alt
import polars as pl

METH_CALLER_TO_NAME = {
    "varlo_1.0": "Varlociraptor α = 1.0",
    "varlo_0.5": "Varlociraptor α = 0.5",
    "varlo_0.1": "Varlociraptor α = 0.1",
    "varlo_0.05": "Varlociraptor α = 0.05",
    "varlo_0.01": "Varlociraptor α = 0.01",
    "varlociraptor": "Varlociraptor",
    "bismark": "Bismark",
    "bsMap": "BSMAPz",
    "methylDackel": "MethylDackel",
    "modkit": "Modkit",
    "pb_CpG_tools": "Pb-CpG-tools",
    "bisSNP": "BisSNP",
}

TOOL_COLORS = {
    "Bismark": "#D81B60",
    "BSMAPz": "#1E88E5",
    "BisSNP": "#FFC107",
    "MethylDackel": "#f0700e",
    "Modkit": "#B42CEA",
    "Pb-CpG-tools": "#8D9279",
    # One teal ramp from dark (α = 1.0) to light (α = 0.01), so that the α levels are distinguishable
    "Varlociraptor α = 1.0": "#00362c",
    "Varlociraptor α = 0.5": "#005345",
    "Varlociraptor α = 0.1": "#00715f",
    "Varlociraptor α = 0.05": "#028d77",
    "Varlociraptor α = 0.01": "#05AA8F",
}

# Same names as in plot_bars_illumina.py
METRIC_TO_NAME = {
    "mape": "Dᵣ",
    "mae": "Dₐ",
    # "binary_discordance": " Dₛ",
}

pl.Config.set_tbl_cols(20)
pl.Config.set_tbl_rows(20)
alt.data_transformers.enable("vegafusion")

sys.stderr = open(snakemake.log[0], "w")


def read_coverage(coverage_file: str, rep_name: str) -> pl.DataFrame:
    return (
        pl.read_csv(
            coverage_file,
            separator="\t",
            has_header=False,
            new_columns=["chromosome", "pos_start", "pos_end", "coverage"],
        )
        .with_columns(
            ((pl.col("pos_end") + pl.col("pos_start")) / 2)
            .cast(pl.Int64)
            .alias("position")
        )
        .select(["chromosome", "position", "coverage"])
        .rename({"coverage": f"coverage_{rep_name}"})
    )


def load_methylation_data(meth_data_path: str, meth_callers: list[str]) -> pl.DataFrame:
    methylation_columns = [
        f"{caller}_methylation_{rep}"
        for caller in meth_callers
        for rep in ["rep1", "rep2"]
    ]
    return (
        pl.scan_parquet(meth_data_path)
        .select(["sample", "chromosome", "position", *methylation_columns])
        .with_columns(pl.col("chromosome").cast(pl.Int64))
        .collect()
    )


"""
Join the coverage of both replicates and compute their minimum per site.
Sites without a coverage entry in either replicate get a null min coverage.
"""
def join_coverage(
    df: pl.DataFrame, coverage_path_rep1: str, coverage_path_rep2: str
) -> pl.DataFrame:
    coverage_rep1 = read_coverage(coverage_path_rep1, "rep1")
    coverage_rep2 = read_coverage(coverage_path_rep2, "rep2")
    return (
        df.join(coverage_rep1, on=["chromosome", "position"], how="left")
        .join(coverage_rep2, on=["chromosome", "position"], how="left")
        .with_columns(
            pl.when(
                pl.col("coverage_rep1").is_not_null()
                & pl.col("coverage_rep2").is_not_null()
            )
            .then(pl.min_horizontal("coverage_rep1", "coverage_rep2"))
            .alias("min_coverage")
        )
    )


def discordance_metrics() -> list[pl.Expr]:
    """Same formulas as compute_replicate_counts in prepare_plot_df.py."""
    rep1 = pl.col("rep1")
    rep2 = pl.col("rep2")
    return [
        pl.len().alias("N"),
        pl.when((rep1 == 0) & (rep2 == 0))
        .then(0.0)
        .otherwise((rep1 - rep2).abs() / pl.max_horizontal(rep1, rep2))
        .mean()
        .alias("mape"),
        (rep1 - rep2).abs().mean().alias("mae"),
        # ((rep1 > 0) != (rep2 > 0)).cast(pl.Float64).mean().alias("binary_discordance"),
    ]


def compute_discordance_vs_loci(
    df: pl.DataFrame, meth_callers: list[str], min_coverages: list[int]
) -> pl.DataFrame:
    metrics = discordance_metrics()
    frames = []
    for caller in meth_callers:
        caller_df = (
            df.select(
                pl.col(f"{caller}_methylation_rep1").alias("rep1"),
                pl.col(f"{caller}_methylation_rep2").alias("rep2"),
                "min_coverage",
            )
            .drop_nulls(["rep1", "rep2"])
        )
        for min_coverage in min_coverages:
            filtered = (
                caller_df
                if min_coverage <= 0
                else caller_df.filter(pl.col("min_coverage") >= min_coverage)
            )
            frames.append(
                filtered.select(metrics).with_columns(
                    pl.lit(caller).alias("caller"),
                    pl.lit(min_coverage).alias("min_coverage"),
                )
            )
    return pl.concat(frames)


def plot_discordance_vs_loci(df: pl.DataFrame):
    color_domain = [
        tool for tool in TOOL_COLORS if tool in df["tool_name"].unique().to_list()
    ]
    color_range = [TOOL_COLORS[tool] for tool in color_domain]

    return (
        alt.Chart(df)
        .mark_line(point=alt.OverlayMarkDef(size=12), strokeWidth=1)
        .encode(
            x=alt.X("N:Q", title="Retained CpG loci (N)"),
            y=alt.Y("value:Q", title="Discordance", scale=alt.Scale(zero=False)),
            order=alt.Order("min_coverage:Q"),
            color=alt.Color(
                "tool_name:N",
                title="Caller",
                scale=alt.Scale(domain=color_domain, range=color_range),
            ),
            tooltip=[
                alt.Tooltip("tool_name:N", title="Caller"),
                alt.Tooltip("min_coverage:Q", title="Min coverage"),
                alt.Tooltip("N:Q", format=","),
                alt.Tooltip("value:Q", title="Discordance", format=".4f"),
            ],
        )
        .properties(width=200, height=200)
        .facet(
            row=alt.Row(
                "metric:N", title=None, sort=list(METRIC_TO_NAME.values())
            ),
            column=alt.Column("sample:N", title=None),
        )
        .resolve_scale(y="independent")
        .configure_header(labelFontSize=14, labelFontWeight="bold")
    )


meth_callers = snakemake.params["meth_callers"]
samples = snakemake.params["samples"]

platform = snakemake.params["platform"]
min_coverages = snakemake.params["min_coverages"]


meth_df = load_methylation_data(snakemake.input["meth_data"], meth_callers)

sample_frames = []
for sample, coverage_path_rep1, coverage_path_rep2 in zip(
    samples, snakemake.input["coverages_rep1"], snakemake.input["coverages_rep2"]
):
    df = join_coverage(
        meth_df.filter(pl.col("sample") == sample), coverage_path_rep1, coverage_path_rep2
    )
    sample_frames.append(
        compute_discordance_vs_loci(df, meth_callers, min_coverages).with_columns(
            # PacBio and Nanopore only have the sample "dummy", so we show the platform instead
            pl.lit(
                platform if sample == "dummy" else sample.replace("_HG002_", "_")
            ).alias("sample")
        )
    )

df_discordance = (
    pl.concat(sample_frames)
    # Make a row per metric instead of a column per metric, so that we can facet the plot by metric. The new column is called "metric" and the values are in a column called "value".
    .unpivot(
        index=["sample", "caller", "min_coverage", "N"],
        on=list(METRIC_TO_NAME.keys()),
        variable_name="metric",
    )
    # Rename methrics and callers to more readable names
    .with_columns(
        pl.col("metric").replace(METRIC_TO_NAME),
        pl.col("caller").replace(METH_CALLER_TO_NAME).alias("tool_name"),
    )
)


chart = plot_discordance_vs_loci(df_discordance)
chart.save(snakemake.output[0], scale_factor=2)
df_discordance.select(
    "sample", "caller", "tool_name", "min_coverage", "N", "metric", "value"
).write_parquet(snakemake.output[1])
