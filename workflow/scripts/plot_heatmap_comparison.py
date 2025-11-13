import pandas as pd
import altair as alt
import sys
import numpy as np
from functools import reduce


sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
alt.data_transformers.enable("vegafusion")


def bin_methylation(series: pd.Series, bin_size: int) -> pd.Series:
    """
    Round methylation values to the nearest bin_size and cast to int.
    """
    return (np.round(series / bin_size) * bin_size).astype(int)


def plot_biases(df):
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

    keep_cols = ["DP", "AF", "SB", "ROB", "RPB", "SCB", "HE", "ALB"]

    # Split rep1
    df_rep1 = df["varlo_format_rep1"].str.split(":", expand=True)
    df_rep1.columns = [f"{col}_rep1" for col in new_cols[: df_rep1.shape[1]]]
    df_rep1 = df_rep1[[f"{col}_rep1" for col in keep_cols]]

    # Split rep2
    df_rep2 = df["varlo_format_rep2"].str.split(":", expand=True)
    df_rep2.columns = [f"{col}_rep2" for col in new_cols[: df_rep2.shape[1]]]
    df_rep2 = df_rep2[[f"{col}_rep2" for col in keep_cols]]

    df_selected = pd.concat(
        [df[["chromosome", "position"]].reset_index(drop=True), df_rep1, df_rep2],
        axis=1,
    )
    # Filter rows without bias info
    bias_cols = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]
    filter_cols = [f"{b}_rep1" for b in bias_cols] + [f"{b}_rep2" for b in bias_cols]
    df_selected = df_selected[df_selected[filter_cols].notna().all(axis=1)]
    df_selected = df_selected[df_selected[filter_cols].ne(".").any(axis=1)]

    # AF-cols to float
    df_selected["AF_rep1"] = df_selected["AF_rep1"].astype(float)
    df_selected["AF_rep2"] = df_selected["AF_rep2"].astype(float)
    df_selected["DP_rep1"] = df_selected["DP_rep1"].astype(int)
    df_selected["DP_rep2"] = df_selected["DP_rep2"].astype(int)

    df_selected["rep1_has_bias"] = (
        df_selected[[f"{b}_rep1" for b in bias_cols]].ne(".").any(axis=1)
    )
    df_selected["rep2_has_bias"] = (
        df_selected[[f"{b}_rep2" for b in bias_cols]].ne(".").any(axis=1)
    )

    # Bias categories
    def categorize(row):
        if row["rep1_has_bias"] and row["rep2_has_bias"]:
            return "Bias both reps"
        elif (row["rep1_has_bias"] and row["AF_rep2"] == 0) or (
            row["rep2_has_bias"] and row["AF_rep1"] == 0
        ):
            return "Bias, AF = 0"
        elif (row["rep1_has_bias"] and row["AF_rep2"] > 0) or (
            row["rep2_has_bias"] and row["AF_rep1"] > 0
        ):
            return "Bias, AF > 0"
        else:
            return "No bias"

    df_selected["category"] = df_selected.apply(categorize, axis=1)

    rep1_long = df_selected.melt(
        id_vars=[
            "chromosome",
            "position",
            "AF_rep1",
            "AF_rep2",
            "category",
            "DP_rep1",
            "DP_rep2",
        ],
        value_vars=[f"{b}_rep1" for b in bias_cols],
        var_name="bias_type",
        value_name="rep1_bias",
    )
    rep1_long["bias_type"] = rep1_long["bias_type"].str.replace("_rep1", "")

    rep2_long = df_selected.melt(
        id_vars=[
            "chromosome",
            "position",
            "AF_rep1",
            "AF_rep2",
            "category",
            "DP_rep1",
            "DP_rep2",
        ],
        value_vars=[f"{b}_rep2" for b in bias_cols],
        var_name="bias_type",
        value_name="rep2_bias",
    )
    rep2_long["bias_type"] = rep2_long["bias_type"].str.replace("_rep2", "")

    df_long = rep1_long.merge(
        rep2_long[["chromosome", "position", "bias_type", "rep2_bias"]],
        on=["chromosome", "position", "bias_type"],
    )

    df_long["rep1_has_bias"] = df_long["rep1_bias"].ne(".")
    df_long["rep2_has_bias"] = df_long["rep2_bias"].ne(".")
    df_long = df_long[df_long["rep1_has_bias"] | df_long["rep2_has_bias"]]

    bias_colors = {
        "SB": "#D81B60",
        "ALB": "#1E88E5",
    }
    # This is only for paper plotting
    if snakemake.output.get("bias_df") is not None:
        df_long.to_parquet(snakemake.output["bias_df"], index=False)
    bias_chart = (
        alt.Chart(df_long)
        .mark_bar()
        .encode(
            x=alt.X("category:N", title="Category"),
            y=alt.Y("count():Q", title="Number of sites"),
            color=alt.Color(
                "bias_type:N",
                title="Bias type",
                scale=alt.Scale(
                    domain=list(bias_colors.keys()),
                    range=list(bias_colors.values()),
                ),
            ),
            tooltip=["category", "bias_type", "count()"],
        )
        .properties(title="Bias Types per Category Across Replicates")
    )
    df_long = df_long[df_long["category"] == "Bias, AF > 0"]
    df_long["AF"] = np.maximum(df_long["AF_rep1"], df_long["AF_rep2"]).round(2)
    df_long = df_long[df_long["DP_rep1"] <= 500]
    df_long = df_long[df_long["DP_rep2"] <= 500]

    af_chart = (
        alt.Chart(df_long)
        .mark_bar()
        .encode(
            x=alt.X("AF:Q", title="Allele frequency (AF)"),
            y=alt.Y("count()"),
        )
        .properties(title="AF of replicates without bias")
    )

    df_long = df_long.melt(
        id_vars=["chromosome", "position"],
        value_vars=["DP_rep1", "DP_rep2"],
        var_name="replicate",
        value_name="DP",
    )

    dp_chart = (
        alt.Chart(df_long)
        .mark_bar()
        .encode(
            x=alt.X("DP:Q", bin=alt.Bin(step=5), title="Depth (DP)"),
            y=alt.Y("count()", title="Anzahl Positionen"),
            color=alt.Color("replicate:N", title="Replicate"),
            # column=alt.Column("replicate:N"),  # nebeneinander (facettiert)
        )
        .properties(title="DP rep1 vs. DP rep2")
    )
    bias_chart = alt.hconcat(bias_chart, af_chart, dp_chart).resolve_scale(
        color="independent"
    )

    return bias_chart, df_long


def plot_count_heatmap(
    df: pd.DataFrame,
    meth_caller: str,
    bin_size: int,
    mapes: dict,
    meth_caller_name: str,
) -> alt.Chart:
    """
    Create a log-scaled heatmap for replicate methylation values.

    Parameters
    ----------
    df : pd.DataFrame
        Binned methylation counts with 'rep1_bin', 'rep2_bin', 'count'.
    meth_caller : str
        Name of the methylation caller.
    bin_size : int
        Bin size for methylation percentages.
    mapes : dict
        Mean absolute percentage errors per methylation caller.

    Returns
    -------
    alt.Chart
        Altair heatmap chart.
    """
    max_count = df["count"].max()
    min_count = 1
    ticks = [min_count, max_count]
    heatmap = (
        alt.Chart(
            df,
            title=alt.Title(
                meth_caller_name,
                subtitle=f"N = {df['count'].sum()} | D = {mapes.get(meth_caller,0)}",
            ),
        )
        .mark_rect()
        .encode(
            x=alt.X(
                "rep1_bin:O",
                sort=list(range(0, 101, bin_size)),
                title="replicate 1",
            ),
            y=alt.Y(
                "rep2_bin:O",
                sort=list(range(100, -1, -bin_size)),
                title="replicate 2",
            ),
            color=alt.Color(
                "count:Q",
                scale=alt.Scale(type="log", scheme="viridis", domain=[1, max_count]),
                legend=alt.Legend(
                    title="Count", orient="right", values=ticks, format=","
                ),
            ),
            tooltip=["rep1_bin", "rep2_bin", "count"],
        )
        .properties(width=200, height=200)
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
            tooltip=["dist_bin", "meth_caller", "rel_count"],
        )
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
            tooltip=["dist", "meth_caller", "cdf"],
        )
    )

    combined_chart = alt.layer(cdf_chart, hist_chart).properties(
        width=700, height=400, title="CDF + Histogram of distances"
    )
    return combined_chart


def compute_replicate_counts(df_dict: dict, bin_size: int, samples):
    """
    Compute replicate counts and binned/unbinned distances for each methylation caller.
    Returns:
        meth_caller_dfs: DataFrame per caller with binned counts
        cdf_dfs: DataFrame per caller for CDFs
        mapes: mean absolute percentage error per caller
    """
    meth_callers = snakemake.params["meth_callers"]
    # Merge across samples if needed
    if isinstance(samples, str):
        samples = [samples]
    df = pd.concat([df_dict[p] for p in samples], ignore_index=True)

    # methyl_cols = [c for c in df.columns if "methylation" in c]
    # df = df.groupby(["chromosome", "position"], as_index=False)[methyl_cols].mean()

    meth_caller_dfs = {}
    cdf_dfs = {}
    mapes = {}

    for meth_caller in meth_callers:
        x_col = f"{meth_caller}_methylation_rep1"
        y_col = f"{meth_caller}_methylation_rep2"

        temp_df = df.dropna(subset=[x_col, y_col])
        if not temp_df.empty:
            mape = (
                np.where(
                    (temp_df[x_col] == 0) & (temp_df[y_col] == 0),
                    0,
                    np.abs(temp_df[x_col] - temp_df[y_col])
                    / temp_df[[x_col, y_col]].max(axis=1),
                ).mean()
                * 100
            )
            mapes[meth_caller] = f"{mape:.2f}%"

        df_temp = temp_df.copy()
        df_temp["rep1_bin"] = bin_methylation(df_temp[x_col], bin_size)
        df_temp["rep2_bin"] = bin_methylation(df_temp[y_col], bin_size)

        # Binned counts for heatmap
        agg_histo = (
            pd.crosstab(df_temp["rep1_bin"], df_temp["rep2_bin"])
            .stack()
            .reset_index(name="count")
        )
        agg_histo["rel_count"] = agg_histo["count"] / agg_histo["count"].sum()
        agg_histo["dist"] = agg_histo.apply(
            lambda row: (
                0
                if max(row["rep1_bin"], row["rep2_bin"]) == 0
                else abs(row["rep1_bin"] - row["rep2_bin"])
                / max(row["rep1_bin"], row["rep2_bin"])
                * 100
            ),
            axis=1,
        )
        agg_histo["dist_bin"] = (agg_histo["dist"] / bin_size).round() * bin_size
        agg_histo["dist_bin"] = agg_histo["dist_bin"].astype(int)
        agg_histo["meth_caller"] = meth_caller
        meth_caller_dfs[meth_caller] = agg_histo

        # Unbinned distances for CDF
        agg_cdf = df_temp.copy()
        # Calculate distance directly using vectorized operations
        max_vals = agg_cdf[[x_col, y_col]].max(axis=1)
        agg_cdf["dist"] = np.where(
            max_vals == 0,
            0,
            (np.abs(agg_cdf[x_col] - agg_cdf[y_col]) / max_vals * 100).round(2),
        )
        agg_cdf = agg_cdf.groupby("dist", as_index=False).agg(count=("dist", "size"))
        agg_cdf["rel_count"] = agg_cdf["count"] / agg_cdf["count"].sum() * 100
        agg_cdf["cdf"] = agg_cdf["count"].cumsum() / agg_cdf["count"].sum()
        agg_cdf["meth_caller"] = meth_caller
        cdf_dfs[meth_caller] = agg_cdf[["dist", "cdf", "meth_caller"]]

    return meth_caller_dfs, cdf_dfs, mapes


# -----------------------------
# Main execution
# -----------------------------
samples = snakemake.params["sample"]
if isinstance(samples, str):
    samples = [samples]

bin_size = snakemake.params["bin_size"]
meth_callers = snakemake.params["meth_callers"]
fdr = snakemake.params["fdr"]
if float(fdr) == 1.0 and snakemake.params.get("paper_plots", False) == True:
    meth_callers = ["varlo"]

# Load HDF5 input
meth_caller_dfs = {}
with pd.HDFStore(snakemake.input[0], mode="r", locking=False) as store:
    for key in store.keys():
        meth_caller_dfs[key.strip("/")] = store[key]

sample_df = pd.concat([meth_caller_dfs[p] for p in samples], ignore_index=True)


############################################## Histogram single Illumina samples
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
                "distance": float(
                    mape_s[m].replace("%", "")
                ),  # Convert "12.3%" -> 12.3
            }
        )
df_summary = pd.DataFrame(results)

colorblind_safe_palette = [
    "#D81B60",
    "#1E88E5",
    "#FFC107",
    "#D35892",
    "#AC3FE6",
]

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
    )
)

labels = (
    alt.Chart(df_summary)
    .mark_text(size=8, dy=-5, color="black")
    .transform_calculate(text_k="datum.number + 'k'")
    .encode(text="text_k:N", x="sample:N", xOffset="meth_caller:N", y="distance:Q")
)

illumina_histo = bars + labels


####################################################### Compute Biases

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

################################################################## Compute heatmaps
replicate_dfs, cdf_dfs, mapes = compute_replicate_counts(
    meth_caller_dfs, bin_size, snakemake.params["sample"]
)


# Create heatmaps for each methylation caller
heatmaps = [
    plot_count_heatmap(replicate_dfs[m], m, bin_size, mapes, meth_caller_to_name[m])
    for m in meth_callers
]
heatmap_df = pd.concat(replicate_dfs.values(), ignore_index=True)
heatmap_plots = alt.hconcat(*heatmaps).resolve_scale(color="independent")

#########################################################

# Save heatmap output
if snakemake.params.get("paper_plots", False) == True:
    df_summary.to_parquet(snakemake.output["histogram_single_samples"])
    bias_df.to_parquet(snakemake.output["bias"])
    heatmap_df.to_parquet(snakemake.output["heatmap"])
else:
    heatmap_plots.save(
        snakemake.output["heatmap"], embed_options={"actions": False}, inline=False
    )
    if snakemake.output.get("histogram_single_samples") is not None:
        illumina_histo.save(
            snakemake.output["histogram_single_samples"],
            embed_options={"actions": False},
            inline=False,
        )
    bias_chart.save(
        snakemake.output["bias"], embed_options={"actions": False}, inline=False
    )
