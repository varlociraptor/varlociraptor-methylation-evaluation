import pandas as pd
import altair as alt
import sys
import numpy as np
import pickle

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
alt.data_transformers.enable("vegafusion")


def bin_methylation(series: pd.Series, bin_size: int) -> pd.Series:
    """Round methylation values to nearest bin_size and cast to int."""

    return (np.round(series / bin_size) * bin_size).astype(int)


def compute_replicate_counts(df_dict: dict, bin_size: int, samples: list):
    meth_callers = snakemake.params["meth_callers"]

    df = pd.concat([df_dict[p] for p in samples], ignore_index=True)

    caller_counts = []
    mape_records = []
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

        mape_records.append({"meth_caller": caller, "mape": mape})

        temp = temp.assign(
            rep1_bin=bin_methylation(temp[rep1], bin_size),
            rep2_bin=bin_methylation(temp[rep2], bin_size),
        )

        counts = (
            pd.crosstab(temp["rep1_bin"], temp["rep2_bin"])
            .stack()
            .reset_index(name="count")
        )

        counts["rel_count"] = counts["count"] / counts["count"].sum()

        max_bin = counts[["rep1_bin", "rep2_bin"]].max(axis=1)
        counts["dist"] = np.where(
            max_bin == 0,
            0,
            np.abs(counts["rep1_bin"] - counts["rep2_bin"]) / max_bin * 100,
        )

        counts["dist_bin"] = (counts["dist"] / bin_size).round().astype(int) * bin_size
        counts["meth_caller"] = caller

        caller_counts.append(counts)

    counts_df = pd.concat(caller_counts, ignore_index=True)
    mapes_df = pd.DataFrame(mape_records)

    return counts_df, mapes_df


samples = snakemake.params["sample"]

if isinstance(samples, str):
    samples = [samples]
plot_type = snakemake.params.get("plot_type")
bin_size = snakemake.params["bin_size"]

meth_callers = snakemake.params["meth_callers"]

# Load HDF5 input
meth_caller_dfs = {}
with pd.HDFStore(snakemake.input[0], mode="r", locking=False) as store:
    for key in store.keys():
        meth_caller_dfs[key.strip("/")] = store[key]
# Combine sample data
sample_df = pd.concat([meth_caller_dfs[p] for p in samples], ignore_index=True)


replicate_dfs, mapes = compute_replicate_counts(meth_caller_dfs, bin_size, samples)


sample_name = snakemake.params["sample_name"].replace("_HG002_", "_")
mapes["sample"] = sample_name
replicate_dfs["sample"] = sample_name

replicate_dfs.to_parquet(snakemake.output["df"], engine="pyarrow")
mapes.to_parquet(snakemake.output["mapes"], engine="pyarrow")
