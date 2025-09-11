import pandas as pd
import altair as alt
import sys
import numpy as np

# Redirect error output to log file
sys.stderr = open(snakemake.log[0], "w")

# Show all columns during debugging
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)


def bin_methylation(series, bin_size=10):
    """Rundet Methylierungswerte auf die nächsten bin_size-Schritte und castet zu int"""
    return (np.round(series / bin_size) * bin_size).astype(int)


def plot_diff_heatmap(method, df_dict):
    df1 = df_dict[method]
    df2 = df_dict["varlo"]
    df_diff = df1.merge(
        df2, on=["rep1_bin", "rep2_bin"], suffixes=(f"_varlo", f"_{method}")
    )
    df_diff["abs_diff"] = (df_diff[f"count_varlo"] - df_diff[f"count_{method}"]).abs()
    df_diff["log2FC"] = np.where(
        (df_diff[f"count_varlo"] != 0) & (df_diff[f"count_{method}"] != 0),
        np.log2(df_diff[f"count_varlo"] / df_diff[f"count_{method}"]).abs(),
        0,
    )

    # Difference Heatmap
    diff_heatmap = (
        alt.Chart(df_diff)
        .mark_rect()
        .encode(
            x=alt.X("rep1_bin:O", sort=list(range(0, 101, 10))),
            y=alt.Y("rep2_bin:O", sort=list(range(100, -1, -10))),
            color=alt.Color(
                "abs_diff:Q",
                scale=alt.Scale(
                    type="log", scheme="viridis", domain=[1, df_diff["abs_diff"].max()]
                ),
            ),
            tooltip=["rep1_bin", "rep2_bin", "abs_diff"],
        )
        .properties(title=f"Abs Difference {method} vs varlo", width=150, height=150)
    )

    # Ratio Heatmap
    ratio_heatmap = (
        alt.Chart(df_diff)
        .mark_rect()
        .encode(
            x=alt.X("rep1_bin:O", sort=list(range(0, 101, 10))),
            y=alt.Y("rep2_bin:O", sort=list(range(100, -1, -10))),
            color=alt.Color(
                "log2FC:Q",
                scale=alt.Scale(scheme="viridis", domain=[0, df_diff["log2FC"].max()]),
            ),
            tooltip=["rep1_bin", "rep2_bin", "log2FC"],
        )
        .properties(title=f"Log2 Fold Change {method}/varlo", width=150, height=150)
    )

    return alt.vconcat(diff_heatmap, ratio_heatmap).resolve_scale(color="independent")


def plot_heatmap_meth_callers(df_dict):
    sample_name = snakemake.params["protocol"]
    df = df_dict[sample_name]

    meth_caller_dfs = {}
    method_plots = []

    # Counts für jede meth_callere berechnen
    for meth_caller in snakemake.params["meth_callers"]:
        print(f"Methode: {meth_caller}", file=sys.stderr)
        x_col = f"{meth_caller}_methylation_rep1"
        y_col = f"{meth_caller}_methylation_rep2"

        df_temp = df.dropna(
            subset=[x_col, y_col]
        ).copy()  # .copy(), um SettingWithCopyWarning zu vermeiden
        df_temp["rep1_bin"] = bin_methylation(df_temp[x_col])
        df_temp["rep2_bin"] = bin_methylation(df_temp[y_col])

        agg = (
            pd.crosstab(df_temp["rep1_bin"], df_temp["rep2_bin"])
            .stack()
            .reset_index(name="count")
        )
        meth_caller_dfs[meth_caller] = agg

        heatmap = (
            alt.Chart(agg)
            .mark_rect()
            .encode(
                x=alt.X(
                    "rep1_bin:O",
                    sort=list(range(0, 101, 10)),
                    title=f"{meth_caller} Rep1",
                ),
                y=alt.Y(
                    "rep2_bin:O",
                    sort=list(range(100, -1, -10)),
                    title=f"{meth_caller} Rep2",
                ),
                color=alt.Color(
                    "count:Q",
                    scale=alt.Scale(
                        type="log", scheme="viridis", domain=[1, agg["count"].max()]
                    ),
                ),
                tooltip=["rep1_bin", "rep2_bin", "count"],
            )
            .properties(title=f"{meth_caller} Heatmap", width=150, height=150)
        )

        method_plots.append(heatmap)
        print("Methode fertig", method_plots, file=sys.stderr)

    return (
        alt.hconcat(*method_plots).resolve_scale(color="independent"),
        meth_caller_dfs,
    )


# Daten aus HDF5 einlesen
meth_caller_dfs = {}
with pd.HDFStore(snakemake.input[0], mode="r", locking=False) as store:
    for key in store.keys():
        meth_caller_dfs[key.strip("/")] = store[key]

# Heatmaps erstellen
chart, meth_caller_dfs = plot_heatmap_meth_callers(meth_caller_dfs)
diff_charts = []
# Differenz-Heatmaps erstellen
for method in snakemake.params["meth_callers"]:
    if method == "varlo":
        continue  # Differenz-Heatmap von varlo zu varlo macht keinen Sinn
    diff_chart = plot_diff_heatmap(method, meth_caller_dfs)
    diff_charts.append(diff_chart)

diff_chart_plots = alt.hconcat(*diff_charts).resolve_scale(color="independent")
chart = alt.vconcat(chart, diff_chart_plots)
# Chart speichern
chart.save(
    snakemake.output[0],
    embed_options={"actions": False},
    inline=False,
)
