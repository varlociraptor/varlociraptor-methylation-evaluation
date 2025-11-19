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


def plot_biases(df):
    """Plot bias categories, AF distribution, and DP distribution."""

    bias_colors = {
        "SB": "#D81B60",
        "ALB": "#1E88E5",
    }
    platform = (
        "Illumina"
        if snakemake.params["platform"] == "Illumina_pe"
        else snakemake.params["platform"]
    )
    bias_label_map = {
        "SB": "Strand Bias",
        "ALB": "Alt Locus Bias",
    }
    df_long = prepare_bias_dataframe(df)
    df_long["bias_type_label"] = df_long["bias_type"].map(bias_label_map)

    bias_chart = (
        alt.Chart(df_long)
        .mark_bar()
        .encode(
            x=alt.X("category:N", title=None),
            y=alt.Y("count():Q", title="Number of sites"),
            color=alt.Color(
                "bias_type_label:N",
                scale=alt.Scale(
                    domain=list(bias_label_map.values()),
                    range=[bias_colors[k] for k in bias_label_map.keys()],
                ),
                #     None if platform != "Nanopore" else alt.Legend(title="Bias type")
                # ),
            ),
            tooltip=["category:N", "count():Q", "bias_type_label:N"],
        )
        .interactive()
        .properties(title=f"{platform} data", height=140)
    )
    # AF chart (AF > 0, DP <= 500)
    df_af = df_long[df_long["category"] == "Bias, AF > 0"]
    df_af["AF"] = np.maximum(df_af["AF_rep1"], df_af["AF_rep2"]).round(2)
    df_af = df_af[(df_af["DP_rep1"] <= 500) & (df_af["DP_rep2"] <= 500)]
    af_chart = (
        alt.Chart(df_af)
        .mark_bar()
        .encode(
            x="AF:Q",
            y="count()",
            color=alt.value("#1E88E5"),
            tooltip=["AF:Q", "count()"],
        )
        .properties(title="AF at loci with bias in other replicates")
        .interactive()
    )

    # DP chart
    df_dp = df_af.melt(
        id_vars=["chromosome", "position"],
        value_vars=["DP_rep1", "DP_rep2"],
        var_name="replicate",
        value_name="DP",
    )
    dp_chart = (
        alt.Chart(df_dp)
        .mark_bar()
        .encode(
            x=alt.X("DP:Q", bin=alt.Bin(step=5), title="Depth (DP)"),
            y=alt.Y(
                "count()",
                title="Coverage per replicate at loci with bias in only one replicate",
            ),
            color=alt.Color(
                "replicate:N",
                scale=alt.Scale(range=list(bias_colors.values())),
            ),
            tooltip=["DP:Q", "count()", "replicate:N"],
        )
        .properties(title="DP rep1 vs. DP rep2")
        .interactive()
    )

    combined = alt.hconcat(bias_chart, af_chart, dp_chart).resolve_scale(
        color="independent"
    )
    return combined, df_long


# -----------------------------
# Main execution
# -----------------------------

samples = snakemake.params["sample"]
if isinstance(samples, str):
    samples = [samples]
plot_type = snakemake.params.get("plot_type")
fdr = snakemake.params["fdr"]


meth_caller_to_name = {
    "varlo": f"Varlociraptor α = {fdr}",
    "bismark": "Bismark",
    "bsMap": "BSMAPz",
    "bisSNP": "BisSNP",
    "methylDackel": "MethylDackel",
    "modkit": "Modkit",
    "pb_CpG_tools": "pb-CpG-tools",
}

# Load HDF5 input
meth_caller_dfs = {}
with pd.HDFStore(snakemake.input[0], mode="r", locking=False) as store:
    for key in store.keys():
        meth_caller_dfs[key.strip("/")] = store[key]
# Combine sample data
sample_df = pd.concat([meth_caller_dfs[p] for p in samples], ignore_index=True)

# -----------------------------
# Compute and plot biases
# -----------------------------
bias_chart, bias_df = plot_biases(
    sample_df[
        [
            "chromosome",
            "position",
            "varlo_format_rep1",
            "varlo_format_rep2",
            "varlo_methylation_rep1",
            "varlo_methylation_rep2",
        ]
    ]
)


if plot_type == "parquet":
    bias_df.to_parquet(snakemake.output[0])
elif plot_type == "pkl":
    with open(snakemake.output[0], "wb") as f:
        pickle.dump(bias_chart, f)
else:
    bias_chart.save(snakemake.output[0], embed_options={"actions": False}, inline=False)
