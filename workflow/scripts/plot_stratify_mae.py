import sys

import altair as alt
import numpy as np
import pandas as pd
import polars as pl

# sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
pl.Config.set_tbl_cols(200)
# pl.Config.set_tbl_rows(200)

alt.data_transformers.enable("vegafusion")

# Am Anfang der Datei statt alt.data_transformers.enable("vegafusion")
# Einfach weglassen oder:
alt.data_transformers.disable_max_rows()

# Dann beim Speichern HTML verwenden statt .save()


def bin_methylation(series: pd.Series, bin_size: int) -> pd.Series:
    """Round methylation values to nearest bin_size and cast to int."""

    return (np.round(series / bin_size) * bin_size).astype(int)


def plot_distances(df):
    distances = pd.DataFrame()

    for caller in snakemake.params["meth_callers"]:
        caller_df = df[
            [f"{caller}_methylation_rep1", f"{caller}_methylation_rep2"]
        ].dropna()
        if caller_df.empty:
            continue
        caller_df["distance"] = (
            np.round(
                (
                    (
                        caller_df[f"{caller}_methylation_rep1"]
                        - caller_df[f"{caller}_methylation_rep2"]
                    )
                    / 5
                )
            )
            * 5
        ).astype(int)
        caller_df["caller"] = caller
        caller_df = caller_df.groupby("distance").size().reset_index(name="count")
        caller_df["absolute"] = caller_df["count"]
        caller_df["relative"] = caller_df["count"] / caller_df["count"].sum()
        caller_df["caller"] = caller
        distances = pd.concat(
            [distances, caller_df[["distance", "caller", "relative", "absolute"]]]
        )

    relative_plot = (
        alt.Chart(
            distances,
            title=alt.Title(
                "relative frequency of distances per caller",
            ),
        )
        .mark_point(size=8, filled=True)
        .encode(
            x=alt.X("distance:Q", title="distance"),
            xOffset="caller:N",
            y=alt.Y("relative:Q", title="relative frequency"),
            color="caller:N",
            tooltip=["caller", "distance", "relative"],
        )
    )

    absolute_plot = (
        alt.Chart(
            distances,
            title=alt.Title(
                "absolute frequency of distances per caller",
            ),
        )
        .mark_point(size=8, filled=True)
        .encode(
            x=alt.X("distance:Q", title="distance"),
            xOffset="caller:N",
            y=alt.Y("absolute:Q", title="absolute frequency"),
            color="caller:N",
            tooltip=["caller", "distance", "absolute"],
        )
    )

    alt.hconcat(relative_plot, absolute_plot).save(snakemake.output.distance_plot)


def compute_replicate_counts(df, bin_size):
    meth_callers = snakemake.params["meth_callers"]
    caller_counts = []
    mape_records = []
    charts = []
    for caller in meth_callers:
        rep1 = f"{caller}_methylation_rep1"
        rep2 = f"{caller}_methylation_rep2"
        temp = df.select(
            "position", rep1, rep2, "coverage_rep1", "coverage_rep2"
        ).drop_nulls()
        # Print rows with highest coverage
        print(temp.sort("coverage_rep1", descending=True).head())
        temp = temp.sample(n=20000, seed=42)
        if temp.is_empty():
            continue
        temp = temp.with_columns((pl.col(rep1) - pl.col(rep2)).abs().alias("mae"))
        print(temp)

        plot_data = temp.to_pandas()

        plot = (
            alt.Chart(plot_data, title=f"Coverage vs MAE for {caller}")
            .mark_point(size=5)
            .encode(
                x=alt.X("coverage_rep1:Q", title="Coverage Rep1"),
                y=alt.Y("coverage_rep2:Q", title="Coverage Rep2"),
                color=alt.Color(
                    "mae:Q", scale=alt.Scale(scheme="viridis"), title="MAE"
                ),
                tooltip=["coverage_rep1", "coverage_rep2", "mae"],
            )
        )
        charts.append(plot)
    charts = alt.hconcat(*charts)
    charts.save(snakemake.output[0], scale_factor=2)


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

# Read and prepare methylation data
df = pl.read_parquet(snakemake.input["meth_data"])
df = (
    df.drop([col for col in df.columns if "format" in col])
    .with_columns(pl.col("chromosome").cast(pl.Int64))
    .filter(pl.col("replicate") == sample)
)

# Read and process coverage data
coverage_rep1 = read_coverage(snakemake.input["coverage_all_01"], "rep1")
coverage_rep2 = read_coverage(snakemake.input["coverage_all_02"], "rep2")

# Join coverage data with methylation data
df = df.join(coverage_rep1, on=["chromosome", "position"], how="left")
df = df.join(coverage_rep2, on=["chromosome", "position"], how="left")

# Filter by selected samples
# df = df.filter(pl.col("replicate").is_in(samples))

print(df)
print(coverage_rep1)
print(coverage_rep2)
# Add column position to covera
# df = df[df["replicate"].is_in(samples)]
# Merge coverage data with methylation data


compute_replicate_counts(df, 5)


# sample_name = snakemake.params["sample_name"].replace("_HG002_", "_")
# distances["sample"] = sample_name
# replicate_dfs["sample"] = sample_name


# replicate_dfs.to_parquet(snakemake.output["df"], engine="pyarrow")
# distances.to_parquet(snakemake.output["distances"], engine="pyarrow")
