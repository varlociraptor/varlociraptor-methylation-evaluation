import sys
from calendar import c
from graphlib import TopologicalSorter
from math import dist
from pickletools import dis

import altair as alt
import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
np.set_printoptions(threshold=np.inf)
alt.data_transformers.enable("vegafusion")


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
    distances = []
    for caller in meth_callers:
        rep1 = f"{caller}_methylation_rep1"
        rep2 = f"{caller}_methylation_rep2"
        temp = df[[rep1, rep2]].dropna()
        if temp.empty:
            continue

        rep1_vals = temp[rep1].to_numpy()
        rep2_vals = temp[rep2].to_numpy()

        denom = np.maximum(rep1_vals, rep2_vals)
        mape = (
            np.where(
                (rep1_vals == 0) & (rep2_vals == 0),
                0,
                np.abs(rep1_vals - rep2_vals) / denom,
            ).mean()
            * 100
        )

        binary_discordance = np.mean((rep1_vals > 0) != (rep2_vals > 0)) * 100

        mae = np.abs(rep1_vals - rep2_vals).mean()

        distances.append(
            {
                "meth_caller": caller,
                "mape": mape,
                "mae": mae,
                "binary_discordance": binary_discordance,
            }
        )

        temp = temp.assign(
            rep1_bin=bin_methylation(temp[rep1], bin_size),
            rep2_bin=bin_methylation(temp[rep2], bin_size),
        )

        counts = (
            pd.crosstab(temp["rep1_bin"], temp["rep2_bin"])
            .stack()
            .reset_index(name="count")
        )

        counts["meth_caller"] = caller

        caller_counts.append(counts)

    counts_df = pd.concat(caller_counts, ignore_index=True)
    distances_df = pd.DataFrame(distances)

    return counts_df, distances_df


samples = snakemake.params["sample"]
plot_type = snakemake.params.get("plot_type")
bin_size = snakemake.params["bin_size"]


df = pd.read_parquet(snakemake.input[0], engine="pyarrow")
df = df[df["sample"].isin(samples)]


plot_distances(df)

replicate_dfs, distances = compute_replicate_counts(df, bin_size)


sample_name = snakemake.params["sample_name"].replace("_HG002_", "_")
distances["sample"] = sample_name
replicate_dfs["sample"] = sample_name

replicate_dfs.to_parquet(snakemake.output["df"], engine="pyarrow")
distances.to_parquet(snakemake.output["distances"], engine="pyarrow")
