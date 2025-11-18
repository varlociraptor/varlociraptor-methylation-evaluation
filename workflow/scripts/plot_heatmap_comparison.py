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
    min_count = 1

    ticks = list(np.logspace(0, np.log10(max_count), num=5).round().astype(int))
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
                legend=alt.Legend(
                    title="Count",
                    orient="right",
                    values=ticks,
                    format=",",
                    tickCount=len(ticks),
                ),
            ),
            tooltip=["rep1_bin:O", "rep2_bin:O", "count:Q"],
        )
        .properties(width=200, height=200)
        .interactive()
    )
    return heatmap


def plot_histogram_cdf(meth_callers, meth_caller_dfs, cdf_dfs):
    """
    Combine histogram of binned distances and CDF of unbinned distances
    for multiple methylation callers.
    """
    # Prepare long-form data
    df_hist_long = pd.concat(meth_caller_dfs.values(), ignore_index=True)
    df_hist_long = df_hist_long.groupby(["dist_bin", "meth_caller"], as_index=False)[
        "rel_count"
    ].sum()
    domain = ["varlo"] + [c for c in meth_callers if c != "varlo"]
    df_cdf_long = pd.concat(cdf_dfs.values(), ignore_index=True)

    # Histogram points
    hist_chart = (
        alt.Chart(df_hist_long)
        .mark_point(size=60, filled=True)
        .encode(
            x=alt.X("dist_bin:Q", title="Discordance"),
            y=alt.Y("rel_count:Q", axis=alt.Axis(format="%")),
            color=alt.Color(
                "meth_caller:N", scale=alt.Scale(scheme="category10", domain=domain)
            ),
            tooltip=["dist_bin:Q", "rel_count:Q", "meth_caller:N"],
        )
        .interactive()
    )

    # CDF line
    cdf_chart = (
        alt.Chart(df_cdf_long)
        .mark_line()
        .encode(
            x=alt.X("dist:Q"),
            y=alt.Y("cdf:Q", title="Fraction of loci", axis=alt.Axis(format="%")),
            color=alt.Color(
                "meth_caller:N", scale=alt.Scale(scheme="category10", domain=domain)
            ),
            tooltip=["dist:Q", "cdf:Q", "meth_caller:N"],
        )
        .interactive()
    )

    combined_chart = alt.layer(cdf_chart, hist_chart).properties(
        width=700, height=400, title="CDF + Histogram of distances"
    )
    return combined_chart


# -----------------------------
# Main execution
# -----------------------------

# Parameters
samples = snakemake.params["sample"]
if isinstance(samples, str):
    samples = [samples]
plot_type = snakemake.params.get("plot_type")
bin_size = snakemake.params["bin_size"]
fdr = snakemake.params["fdr"]
meth_callers = (
    ["varlo"]
    if float(fdr) == 0.01 and plot_type in {"parquet", "pkl"}
    else snakemake.params["meth_callers"]
)


if float(fdr) == 0.01 and snakemake.params.get("paper_plots", False) == True:
    meth_callers = ["varlo"]

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
    if plot_type == "parquet":
        bias_df.to_parquet(snakemake.output["bias"])
    elif plot_type == "pkl":
        with open(snakemake.output["bias"], "wb") as f:
            pickle.dump(bias_chart, f)
    else:
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
    if plot_type == "parquet":
        pd.DataFrame().to_parquet(snakemake.output["heatmap"])
    elif plot_type == "pkl":
        with open(snakemake.output["heatmap"], "wb") as f:
            pickle.dump(heatmap_plots, f)
    else:
        heatmap_plots.save(
            snakemake.output["heatmap"], embed_options={"actions": False}, inline=False
        )


# -----------------------------
# Summary histogram for single samples
# -----------------------------
results = []  # List for collecting rows
meth_caller_to_name = {
    "varlo": f"Varlociraptor α = {fdr}",
    "bismark": "Bismark",
    "bsMap": "BSMAPz",
    "methylDackel": "MethylDackel",
    "modkit": "Modkit",
    "pb_CpG_tools": "pb-CpG-tools",
}
for s in samples:
    for m in meth_callers:
        rep_df_s, cdf_s, mape_s = compute_replicate_counts(
            meth_caller_dfs, bin_size, [s]
        )

        results.append(
            {
                "sample": s,
                "meth_caller": meth_caller_to_name[m],
                "number": int(str(rep_df_s[m]["count"].sum())[:3]),
                "distance": float(mape_s[m]),
            }
        )
df_summary = pd.DataFrame(results)

colorblind_safe_palette = [
    "#D81B60",
    "#1E88E5",
    "#FFC107",
    "#05AA8F",
    "#004D40",
]
df_summary["sample"] = df_summary["sample"].str.replace("_HG002_", "_", regex=False)

bars = (
    alt.Chart(df_summary)
    .mark_bar()
    .encode(
        x=alt.X("sample:N"),
        xOffset="meth_caller:N",
        y=alt.Y("distance:Q", title="Discordance"),
        color=alt.Color(
            "meth_caller:N",
            title="Methylation caller",
            scale=alt.Scale(range=colorblind_safe_palette),
        ),
        tooltip=["sample:N", "meth_caller:N", "distance:Q", "number:Q"],
    )
    .interactive()
)

labels = (
    alt.Chart(df_summary)
    .mark_text(size=8, dy=-5, color="black")
    .transform_calculate(text_k="datum.number + 'k'")
    .encode(
        text="text_k:N",
        x="sample:N",
        xOffset="meth_caller:N",
        y="distance:Q",
        tooltip=["sample:N", "meth_caller:N", "distance:Q", "number:Q"],
    )
    .interactive()
)

illumina_histo = bars + labels

if snakemake.output.get("bar_plot_single_samples") is not None:
    if plot_type == "parquet":
        df_summary.to_parquet(snakemake.output["bar_plot_single_samples"])
    elif plot_type == "pkl":
        with open(snakemake.output["bar_plot_single_samples"], "wb") as f:
            pickle.dump(illumina_histo, f)
    else:
        illumina_histo.save(
            snakemake.output["bar_plot_single_samples"],
            embed_options={"actions": False},
            inline=False,
        )
