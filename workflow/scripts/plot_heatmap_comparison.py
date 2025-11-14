import pandas as pd
import altair as alt
import sys
import numpy as np

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
alt.data_transformers.enable("vegafusion")


def bin_methylation(series: pd.Series, bin_size: int) -> pd.Series:
    """Round methylation values to nearest bin_size and cast to int."""
    return (np.round(series / bin_size) * bin_size).astype(int)


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


def plot_biases(df: pd.DataFrame) -> tuple:
    """Generate bias charts and prepare dataframe for plotting."""
    df_long = prepare_bias_dataframe(df)
    colorblind_safe_palette = ["#D81B60", "#1E88E5", "#FFC107", "#D35892", "#AC3FE6"]
    # Bias count chart
    bias_chart = (
        alt.Chart(df_long)
        .mark_bar()
        .encode(
            x="category:N",
            y="count():Q",
            color=alt.Color(
                "bias_type:N",
                scale=alt.Scale(range=colorblind_safe_palette),
            ),
            tooltip=["category", "bias_type", "count()"],
        )
        .properties(title="Bias Types per Category Across Replicates")
    )

    # AF chart (AF > 0, DP <= 500)
    df_af = df_long[df_long["category"] == "Bias, AF > 0"]
    df_af["AF"] = np.maximum(df_af["AF_rep1"], df_af["AF_rep2"]).round(2)
    df_af = df_af[(df_af["DP_rep1"] <= 500) & (df_af["DP_rep2"] <= 500)]
    af_chart = (
        alt.Chart(df_af)
        .mark_bar()
        .encode(x="AF:Q", y="count()")
        .properties(title="AF at loci with bias in other replicates")
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
                scale=alt.Scale(range=colorblind_safe_palette),
            ),
        )
        .properties(title="DP rep1 vs. DP rep2")
    )

    combined = alt.hconcat(bias_chart, af_chart, dp_chart).resolve_scale(
        color="independent"
    )
    return combined, df_long


def compute_replicate_counts(df_dict: dict, bin_size: int, samples: list):
    """Compute binned counts, distances, and CDFs per methylation caller."""
    meth_callers = snakemake.params["meth_callers"]
    df = pd.concat([df_dict[p] for p in samples], ignore_index=True)

    meth_caller_dfs, cdf_dfs, mapes = {}, {}, {}
    for caller in meth_callers:
        rep1_col = f"{caller}_methylation_rep1"
        rep2_col = f"{caller}_methylation_rep2"
        temp = df.dropna(subset=[rep1_col, rep2_col])
        if temp.empty:
            continue

        # MAPE
        mapes[caller] = (
            np.where(
                (temp[rep1_col] == 0) & (temp[rep2_col] == 0),
                0,
                np.abs(temp[rep1_col] - temp[rep2_col])
                / temp[[rep1_col, rep2_col]].max(axis=1),
            ).mean()
            * 100
        )

        temp["rep1_bin"] = bin_methylation(temp[rep1_col], bin_size)
        temp["rep2_bin"] = bin_methylation(temp[rep2_col], bin_size)

        # Binned counts
        counts = (
            pd.crosstab(temp["rep1_bin"], temp["rep2_bin"])
            .stack()
            .reset_index(name="count")
        )
        counts["rel_count"] = counts["count"] / counts["count"].sum()
        counts["dist"] = np.where(
            counts[["rep1_bin", "rep2_bin"]].max(axis=1) == 0,
            0,
            np.abs(counts["rep1_bin"] - counts["rep2_bin"])
            / counts[["rep1_bin", "rep2_bin"]].max(axis=1)
            * 100,
        )
        counts["dist_bin"] = (counts["dist"] / bin_size).round() * bin_size
        counts["dist_bin"] = counts["dist_bin"].astype(int)
        counts["meth_caller"] = caller
        meth_caller_dfs[caller] = counts

        # Unbinned distances (CDF)
        cdf = temp.copy()
        max_vals = cdf[[rep1_col, rep2_col]].max(axis=1)
        cdf["dist"] = np.where(
            max_vals == 0, 0, np.abs(cdf[rep1_col] - cdf[rep2_col]) / max_vals * 100
        )
        cdf = cdf.groupby("dist", as_index=False).agg(count=("dist", "size"))
        cdf["rel_count"] = cdf["count"] / cdf["count"].sum() * 100
        cdf["cdf"] = cdf["count"].cumsum() / cdf["count"].sum()
        cdf["meth_caller"] = caller
        cdf_dfs[caller] = cdf[["dist", "cdf", "meth_caller"]]

    return meth_caller_dfs, cdf_dfs, mapes


def plot_count_heatmap(
    df: pd.DataFrame,
    meth_caller: str,
    bin_size: int,
    mapes: dict,
    meth_caller_name: str,
) -> alt.Chart:
    """Log-scaled heatmap for replicate methylation counts."""
    max_count = df["count"].max()
    heatmap = (
        alt.Chart(
            df,
            title=alt.Title(
                meth_caller_name,
                subtitle=f"N = {df['count'].sum()} | D = {mapes.get(meth_caller,0):.2f}%",
            ),
        )
        .mark_rect()
        .encode(
            x=alt.X(
                "rep1_bin:O", sort=list(range(0, 101, bin_size)), title="Replicate 1"
            ),
            y=alt.Y(
                "rep2_bin:O", sort=list(range(100, -1, -bin_size)), title="Replicate 2"
            ),
            color=alt.Color(
                "count:Q",
                scale=alt.Scale(type="log", scheme="viridis", domain=[1, max_count]),
                legend=alt.Legend(title="Count", orient="right"),
            ),
            tooltip=["rep1_bin", "rep2_bin", "count"],
        )
        .properties(width=200, height=200)
    )
    return heatmap


# -----------------------------
# Main execution
# -----------------------------

# Parameters
samples = snakemake.params["sample"]
if isinstance(samples, str):
    samples = [samples]

bin_size = snakemake.params["bin_size"]
meth_callers = snakemake.params["meth_callers"]
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

if snakemake.output.get("bias") is not None:
    bias_chart.save(
        snakemake.output["bias"], embed_options={"actions": False}, inline=False
    )

# -----------------------------
# Compute replicate counts / heatmaps
# -----------------------------
replicate_dfs, cdf_dfs, mapes = compute_replicate_counts(
    meth_caller_dfs, bin_size, samples
)

heatmaps = [
    plot_count_heatmap(replicate_dfs[m], m, bin_size, mapes, meth_caller_to_name[m])
    for m in meth_callers
]
heatmap_plots = alt.hconcat(*heatmaps).resolve_scale(color="independent")

if snakemake.output.get("heatmap") is not None:
    heatmap_plots.save(
        snakemake.output["heatmap"], embed_options={"actions": False}, inline=False
    )

# -----------------------------
# Summary histogram for single samples
# -----------------------------
results = []
for s in samples:
    for m in meth_callers:
        rep_df_s, cdf_s, mape_s = compute_replicate_counts(
            meth_caller_dfs, bin_size, [s]
        )
        total_count = int(str(rep_df_s[m]["count"].sum())[:3])
        results.append(
            {
                "sample": s,
                "meth_caller": meth_caller_to_name[m],
                "number": total_count,
                "distance": (
                    float(mape_s[m])
                    if isinstance(mape_s[m], (int, float))
                    else float(str(mape_s[m]).replace("%", ""))
                ),
            }
        )

df_summary = pd.DataFrame(results)

colorblind_safe_palette = ["#D81B60", "#1E88E5", "#FFC107", "#D35892", "#AC3FE6"]

bars = (
    alt.Chart(df_summary)
    .mark_bar()
    .encode(
        x=alt.X("sample:N", axis=alt.Axis(labelAngle=-30)),
        xOffset="meth_caller:N",
        y=alt.Y("distance:Q", title="Discordance"),
        color=alt.Color(
            "meth_caller:N",
            scale=alt.Scale(range=colorblind_safe_palette),
            title="Methylation caller",
        ),
    )
)

labels = (
    alt.Chart(df_summary)
    .mark_text(size=8, dy=-5, color="black")
    .transform_calculate(text_k="datum.number + 'k'")
    .encode(text="text_k:N", x="sample:N", xOffset="meth_caller:N", y="distance:Q")
)

illumina_histo = bars + labels

if snakemake.output.get("bar_plot_single_samples") is not None:
    illumina_histo.save(
        snakemake.output["bar_plot_single_samples"],
        embed_options={"actions": False},
        inline=False,
    )
