from functools import total_ordering
from mailbox import linesep
from pdb import line_prefix

import altair as alt
import pandas as pd
import polars as pl

# sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
pl.Config.set_tbl_cols(200)
pl.Config.set_tbl_rows(20)

sys.stderr = open(snakemake.log[0], "w")

alt.data_transformers.enable("vegafusion")
alt.data_transformers.disable_max_rows()


def bin_coverages(expr, bin_size):
    return ((expr / bin_size).round(0) * bin_size).cast(pl.Int64)


def plot_cov(df, caller, ticks_axis, feature):
    x = (
        (df[feature] * df["count"]).sum() / df["count"].sum()
        if feature != "count"
        else df["count"].sum()
    )
    plot = (
        alt.Chart(df, title=f"{caller} - {feature}: {x:.2f}")
        .mark_rect()
        .encode(
            x=alt.X(
                "rep1_bin:O", axis=alt.Axis(values=ticks_axis), title="Replicate 1"
            ),
            y=alt.Y(
                "rep2_bin:O",
                axis=alt.Axis(values=ticks_axis),
                sort="descending",
                title="Replicate 2",
            ),
            color=alt.Color(
                f"{feature}:Q",
            ),
            tooltip=["rep1_bin", "rep2_bin", "mean_mae", "mean_mape", "count"],
        )
        .properties(width=200, height=200)
    )
    return plot


def stratify_cov(df):
    mae_plots, mape_plots, count_plots = [], [], []
    for caller in meth_callers:
        df_temp = df.filter(pl.col("tool") == caller).sort("rep1_bin", "rep2_bin")
        rep1_cutoff = df_temp.select(pl.col("rep1_bin").quantile(quantile)).item()
        rep2_cutoff = df_temp.select(pl.col("rep2_bin").quantile(quantile)).item()
        df_temp = (
            df_temp.filter(
                (pl.col("rep1_bin") <= rep1_cutoff)
                & (pl.col("rep2_bin") <= rep2_cutoff)
            )
            .group_by(["rep1_bin", "rep2_bin"])
            .agg(
                pl.col("mae").mean().alias("mean_mae"),
                pl.col("mape").mean().alias("mean_mape"),
                pl.len().alias("count"),
            )
            .to_pandas()
        )
        ticks_axis = list(range(0, int(max(rep1_cutoff, rep2_cutoff)) + 5, 5))
        mae_plots.append(plot_cov(df_temp, caller, ticks_axis, "mean_mae"))
        mape_plots.append(plot_cov(df_temp, caller, ticks_axis, "mean_mape"))
        count_plots.append(plot_cov(df_temp, caller, ticks_axis, "count"))
    mae_plot = alt.hconcat(*mae_plots)
    mape_plot = alt.hconcat(*mape_plots)
    count_plot = alt.hconcat(*count_plots)
    chart = alt.vconcat(mae_plot, mape_plot, count_plot).resolve_scale(
        color="independent"
    )
    return chart


def read_coverage(coverage_file: str, rep_name: str) -> pl.DataFrame:
    """Read and process coverage file, renaming coverage column appropriately."""
    return (
        pl.read_csv(
            coverage_file,
            separator="\t",
            has_header=False,
            new_columns=["chromosome", "pos_start", "pos_end", "coverage"],
        )
        .with_columns(
            ((pl.col("pos_end") + pl.col("pos_start")) / 2)
            .alias("position")
            .cast(pl.Int64)
        )
        .select(["chromosome", "position", "coverage"])
        .rename({"coverage": f"coverage_{rep_name}"})
    )


# Main execution
sample = snakemake.params["sample"]
plot_type = snakemake.params.get("plot_type")
bin_size = snakemake.params.get("bin_size", 1)
quantile = snakemake.params.get("quantile", 1.0)
meth_callers = snakemake.params["meth_callers"]


# Read and prepare methylation data
df = pl.read_parquet(snakemake.input["meth_data"])
df = (
    df.drop([col for col in df.columns if "format" in col])
    .with_columns(pl.col("chromosome").cast(pl.Int64))
    .filter(pl.col("replicate") == sample)
)

# Read and process coverage data
coverage_rep1 = read_coverage(snakemake.input["coverage_01"], "rep1")
coverage_rep2 = read_coverage(snakemake.input["coverage_02"], "rep2")


# Join coverage data with methylation data
df_cov_all = df.join(coverage_rep1, on=["chromosome", "position"], how="left")
df_cov_all = df_cov_all.join(coverage_rep2, on=["chromosome", "position"], how="left")
df_cov_all = df_cov_all.with_columns(
    [
        bin_coverages(pl.col("coverage_rep1"), bin_size).alias("rep1_bin"),
        bin_coverages(pl.col("coverage_rep2"), bin_size).alias("rep2_bin"),
    ]
)
# Compute MAPE and MAE
df_cov_all = df_cov_all.with_columns(
    [
        expr
        for caller in meth_callers
        for expr in (
            (
                pl.when(
                    pl.max_horizontal(
                        pl.col(f"{caller}_methylation_rep1"),
                        pl.col(f"{caller}_methylation_rep2"),
                    )
                    != 0
                )
                .then(
                    (
                        (
                            pl.col(f"{caller}_methylation_rep1")
                            - pl.col(f"{caller}_methylation_rep2")
                        ).abs()
                        * 100
                    )
                    / pl.max_horizontal(
                        pl.col(f"{caller}_methylation_rep1"),
                        pl.col(f"{caller}_methylation_rep2"),
                    )
                )
                .otherwise(0)
                .alias(f"{caller}_mape")
            ),
            (
                (
                    pl.col(f"{caller}_methylation_rep1")
                    - pl.col(f"{caller}_methylation_rep2")
                )
                .abs()
                .alias(f"{caller}_mae")
            ),
        )
    ]
)
dfs = []

for caller in meth_callers:
    df_tmp = df_cov_all.select(
        [
            "chromosome",
            "position",
            pl.col(f"{caller}_methylation_rep1").alias("methylation_rep1"),
            pl.col(f"{caller}_methylation_rep2").alias("methylation_rep2"),
            pl.col(f"{caller}_mae").alias("mae"),
            pl.col(f"{caller}_mape").alias("mape"),
            pl.col("rep1_bin").alias("rep1_bin"),
            pl.col("rep2_bin").alias("rep2_bin"),
        ]
    ).with_columns(
        pl.lit(caller).alias("tool"),
        bin_coverages(
            pl.max_horizontal(
                pl.col("methylation_rep1"),
                pl.col("methylation_rep2"),
            ),
            5,
        ).alias("meth_bin"),
    )
    dfs.append(df_tmp)

df_long = pl.concat(dfs)
# Plot df_long as a line chart with color = tool, y-axis = mae and x-axis = meth_bin with altair
df_meth = df_long.group_by("tool", "meth_bin").agg(
    pl.col("mae").mean().alias("mae"),
    pl.col("mape").mean().alias("mape"),
    pl.len().alias("count"),
)
#
mae_by_meth_bin = (
    alt.Chart(df_meth)
    .mark_line()
    .encode(
        x="meth_bin",
        y="mae",
        color="tool",
    )
)
mape_by_meth_bin = (
    alt.Chart(df_meth)
    .mark_line()
    .encode(
        x="meth_bin",
        y="mape",
        color="tool",
    )
)

count_by_meth_bin = (
    alt.Chart(df_meth)
    .mark_line()
    .encode(
        x="meth_bin",
        y="count",
        color="tool",
    )
)

line_chart = alt.vconcat(mae_by_meth_bin, mape_by_meth_bin, count_by_meth_bin)


strat_by_cov = stratify_cov(df_long)
chart = alt.hconcat(strat_by_cov, line_chart)
chart.save(snakemake.output[0], scale_factor=2)
