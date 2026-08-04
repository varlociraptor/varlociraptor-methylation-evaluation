import sys

import altair as alt
import pandas as pd
import polars as pl

METH_CALLER_TO_NAME = {
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
    "Varlociraptor α = 0.1": "#004D40",
    "Varlociraptor α = 0.05": "#126e5f",
    "Varlociraptor α = 0.01": "#05AA8F",
}

COVERAGE_COLUMNS = ["chromosome", "pos_start", "pos_end", "coverage"]
MIN_COVERAGE_BIN_SIZE = 1


pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)
pd.set_option("display.expand_frame_repr", False)
pd.set_option("display.max_rows", None)
pl.Config.set_tbl_cols(20)
pl.Config.set_tbl_rows(20)
alt.data_transformers.enable("vegafusion")


sys.stderr = open(snakemake.log[0], "w")


def bin_values(expr: pl.Expr, bin_size: int) -> pl.Expr:
    return ((expr / bin_size).round(0) * bin_size).cast(pl.Int64)


def read_coverage(coverage_file: str, rep_name: str) -> pl.DataFrame:
    """Read a BED-like coverage file and return per-position coverage for one replicate."""
    return (
        pl.read_csv(
            coverage_file,
            separator="\t",
            has_header=False,
            new_columns=COVERAGE_COLUMNS,
        )
        .with_columns(
            ((pl.col("pos_end") + pl.col("pos_start")) / 2)
            .cast(pl.Int64)
            .alias("position")
        )
        .select(["chromosome", "position", "coverage"])
        .rename({"coverage": f"coverage_{rep_name}"})
    )


def load_methylation_data(meth_data_path: str, sample: str) -> pl.DataFrame:
    """Load methylation calls, drop unused format columns, and filter to a sample if needed."""
    df = pl.read_parquet(meth_data_path).with_columns(pl.col("chromosome").cast(pl.Int64))
    df = df.drop([col for col in df.columns if "format" in col])
    if sample != "all_samples":
        df = df.filter(pl.col("sample") == sample)
    return df


def join_coverage(
    df: pl.DataFrame, coverage_path_rep1: str, coverage_path_rep2: str
) -> pl.DataFrame:
    coverage_rep1 = read_coverage(coverage_path_rep1, "rep1")
    coverage_rep2 = read_coverage(coverage_path_rep2, "rep2")
    return df.join(coverage_rep1, on=["chromosome", "position"], how="left").join(
        coverage_rep2, on=["chromosome", "position"], how="left"
    )


def build_long_format(
    df_cov_all: pl.DataFrame, meth_callers: list[str], coverage_bin_size: int
) -> pl.DataFrame:
    """Reshape wide per-caller columns into a long, per-tool table with binned coverage."""
    per_tool_frames = [
        df_cov_all.select(
            "chromosome",
            "position",
            pl.col(f"{caller}_methylation_rep1").alias("methylation_rep1"),
            pl.col(f"{caller}_methylation_rep2").alias("methylation_rep2"),
            "coverage_rep1",
            "coverage_rep2",
        ).with_columns(pl.lit(caller).alias("tool"))
        for caller in meth_callers
    ]

    return (
        pl.concat(per_tool_frames)
        .drop_nulls()
        .with_columns(
            pl.min_horizontal("coverage_rep1", "coverage_rep2").alias("min_coverage"),
        )
    )


def compute_min_coverage_distribution(df_long: pl.DataFrame) -> pl.DataFrame:
    """Bin the min coverage per tool and compute the fraction of sites in each bin."""
    df_min_coverage = (
        df_long.with_columns(
            bin_values(pl.col("min_coverage"), MIN_COVERAGE_BIN_SIZE).alias(
                "min_coverage_bin"
            )
        )
        .group_by("tool", "min_coverage_bin")
        .agg(pl.len().alias("count"))
        .sort("tool", "min_coverage_bin")
        .with_columns(
            (pl.col("count") / pl.sum("count").over("tool")).alias("fraction")
        )
        .with_columns(
            pl.col("fraction").cum_sum().over("tool").alias("cumulative_fraction")
        )
    )
    return df_min_coverage.with_columns(
        pl.col("tool").replace(METH_CALLER_TO_NAME).alias("tool_name")
    )



def plot_min_coverage_vs_fraction(
    df_min_coverage: pl.DataFrame
) -> alt.Chart:
    color_domain = sorted(df_min_coverage["tool_name"].unique())
    color_range = [TOOL_COLORS[tool] for tool in color_domain]
    df_min_coverage = df_min_coverage.filter(pl.col("min_coverage_bin") <= 80)

    return (
        alt.Chart(df_min_coverage)
        .mark_line()
        .encode(
            x=alt.X("min_coverage_bin:Q", title="Min Coverage"),
            y=alt.Y("cumulative_fraction:Q", title="Fraction of called CpG loci (N)"),
            color=alt.Color(
                "tool_name",
                title="Caller",
                scale=alt.Scale(
                    domain=color_domain,
                    range=color_range,
                ),
            ),
            strokeWidth=alt.value(1),
        )
    ).properties(width=200, height=200).configure_header(labelFontSize=14, labelFontWeight="bold")




sample = snakemake.params["sample"]
coverage_bin_size = snakemake.params.get("bin_size", 1)
meth_callers = snakemake.params["meth_callers"]

df = load_methylation_data(snakemake.input["meth_data"], sample)
df_cov_all = join_coverage(
    df, snakemake.input["coverage_01"], snakemake.input["coverage_02"]
)
df_long = build_long_format(df_cov_all, meth_callers, coverage_bin_size)
df_min_coverage = compute_min_coverage_distribution(df_long)

chart = plot_min_coverage_vs_fraction(df_min_coverage)
chart.save(snakemake.output[0], scale_factor=2)
df_min_coverage.write_parquet(snakemake.output[1])
