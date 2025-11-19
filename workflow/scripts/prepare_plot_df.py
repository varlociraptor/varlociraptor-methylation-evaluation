import pandas as pd
import altair as alt
import sys
import numpy as np
import pickle

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
alt.data_transformers.enable("vegafusion")


def split_varlo_format(df: pd.DataFrame, rep: str, keep_cols=None) -> pd.DataFrame:
    """
    Split the varlo_format string into separate columns and select relevant ones.
    """
    if keep_cols is None:
        keep_cols = ["DP", "AF", "SB", "ROB", "RPB", "SCB", "HE", "ALB"]
    new_cols = [
        "DP",
        "AF",
        "SAOBS",
        "SROBS",
        "OBS",
        "OOBS",
        "SB",
        "ROB",
        "RPB",
        "SCB",
        "HE",
        "ALB",
        "AFD",
    ]

    df_split = df[f"varlo_format_{rep}"].str.split(":", expand=True)
    df_split.columns = [f"{col}_{rep}" for col in new_cols[: df_split.shape[1]]]
    return df_split[[f"{col}_{rep}" for col in keep_cols]]


def categorize_bias(row):
    """Assign bias category based on replicates and allele frequency."""
    rep1_has_bias, rep2_has_bias = row["rep1_has_bias"], row["rep2_has_bias"]
    if rep1_has_bias and rep2_has_bias:
        return "Bias both reps"
    elif (rep1_has_bias and row["AF_rep2"] == 0) or (
        rep2_has_bias and row["AF_rep1"] == 0
    ):
        return "Bias, AF = 0"
    elif (rep1_has_bias and row["AF_rep2"] > 0) or (
        rep2_has_bias and row["AF_rep1"] > 0
    ):
        return "Bias, AF > 0"
    return "No bias"


def prepare_bias_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Process Varlociraptor data to identify biases and reshape for plotting."""
    df_rep1 = split_varlo_format(df, "rep1")
    df_rep2 = split_varlo_format(df, "rep2")
    df_selected = pd.concat(
        [df[["chromosome", "position"]].reset_index(drop=True), df_rep1, df_rep2],
        axis=1,
    )
    bias_cols = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]
    filter_cols = [f"{b}_rep1" for b in bias_cols] + [f"{b}_rep2" for b in bias_cols]
    df_selected = df_selected[df_selected[filter_cols].notna().all(axis=1)]
    df_selected = df_selected[df_selected[filter_cols].ne(".").any(axis=1)]

    df_selected[["AF_rep1", "AF_rep2"]] = df_selected[["AF_rep1", "AF_rep2"]].astype(
        float
    )
    df_selected[["DP_rep1", "DP_rep2"]] = df_selected[["DP_rep1", "DP_rep2"]].astype(
        int
    )

    df_selected["rep1_has_bias"] = (
        df_selected[[f"{b}_rep1" for b in bias_cols]].ne(".").any(axis=1)
    )
    df_selected["rep2_has_bias"] = (
        df_selected[[f"{b}_rep2" for b in bias_cols]].ne(".").any(axis=1)
    )

    df_selected["category"] = df_selected.apply(categorize_bias, axis=1)

    # Reshape for plotting
    long_dfs = []
    for rep in ["rep1", "rep2"]:
        temp = df_selected.melt(
            id_vars=[
                "chromosome",
                "position",
                "AF_rep1",
                "AF_rep2",
                "category",
                "DP_rep1",
                "DP_rep2",
            ],
            value_vars=[f"{b}_{rep}" for b in bias_cols],
            var_name="bias_type",
            value_name=f"{rep}_bias",
        )
        temp["bias_type"] = temp["bias_type"].str.replace(f"_{rep}", "")
        long_dfs.append(temp)

    df_long = long_dfs[0].merge(
        long_dfs[1][["chromosome", "position", "bias_type", "rep2_bias"]],
        on=["chromosome", "position", "bias_type"],
    )

    df_long = df_long[df_long["rep1_bias"].ne(".") | df_long["rep2_bias"].ne(".")]

    return df_long


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

        # --- MAPE -----------------------------------------------------------
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

        # --- Binning --------------------------------------------------------
        temp = temp.assign(
            rep1_bin=bin_methylation(temp[rep1], bin_size),
            rep2_bin=bin_methylation(temp[rep2], bin_size),
        )

        # --- Crosstab -------------------------------------------------------
        counts = (
            pd.crosstab(temp["rep1_bin"], temp["rep2_bin"])
            .stack()
            .reset_index(name="count")
        )

        counts["rel_count"] = counts["count"] / counts["count"].sum()

        # Normalized distance %
        max_bin = counts[["rep1_bin", "rep2_bin"]].max(axis=1)
        counts["dist"] = np.where(
            max_bin == 0,
            0,
            np.abs(counts["rep1_bin"] - counts["rep2_bin"]) / max_bin * 100,
        )

        # Distance binning
        counts["dist_bin"] = (counts["dist"] / bin_size).round().astype(int) * bin_size
        counts["meth_caller"] = caller

        caller_counts.append(counts)

    # Combine all callers into single DFs
    counts_df = pd.concat(caller_counts, ignore_index=True)
    mapes_df = pd.DataFrame(mape_records)

    return counts_df, mapes_df


# Parameters
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


# -----------------------------
# Compute replicate counts / heatmaps
# -----------------------------
replicate_dfs, mapes = compute_replicate_counts(meth_caller_dfs, bin_size, samples)
# combined_counts_df = pd.concat(replicate_dfs.values(), ignore_index=True)

sample_name = snakemake.params["sample_name"].replace("_HG002_", "_")
mapes["sample"] = sample_name
replicate_dfs["sample"] = sample_name

replicate_dfs.to_parquet(snakemake.output["df"], engine="pyarrow")
mapes.to_parquet(snakemake.output["mapes"], engine="pyarrow")
