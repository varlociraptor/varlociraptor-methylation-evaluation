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


def plot_heatmap_meth_callers(df_dict):
    sample_name = snakemake.params["protocol"]
    df = df_dict[sample_name]

    meth_caller_dfs = {}

    # Counts für jede meth_callere berechnen
    for meth_caller in snakemake.params["meth_callers"]:
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

    # Wir gehen davon aus: genau 2 meth_calleren
    m1, m2 = snakemake.params["meth_callers"]
    df_diff = meth_caller_dfs[m1].merge(
        meth_caller_dfs[m2], on=["rep1_bin", "rep2_bin"], suffixes=(f"_{m1}", f"_{m2}")
    )
    df_diff["abs_diff"] = (df_diff[f"count_{m1}"] - df_diff[f"count_{m2}"]).abs()
    df_diff["log2FC"] = np.where(
        (df_diff[f"count_{m1}"] != 0) & (df_diff[f"count_{m2}"] != 0),
        np.log2(df_diff[f"count_{m1}"] / df_diff[f"count_{m2}"]).abs(),
        0,
    )
    # Original Heatmaps
    charts = []
    for meth_caller, agg in meth_caller_dfs.items():
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
        charts.append(heatmap)

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
        .properties(title=f"Abs Difference {m1} vs {m2}", width=150, height=150)
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
        .properties(title=f"Log2 Fold Change {m1}/{m2}", width=150, height=150)
    )

    tool_charts = alt.hconcat(*charts).resolve_scale(color="independent")
    diff_charts = alt.hconcat(diff_heatmap, ratio_heatmap).resolve_scale(
        color="independent"
    )
    return alt.vconcat(tool_charts, diff_charts).resolve_scale(color="independent")


# Daten aus HDF5 einlesen
meth_caller_dfs = {}
with pd.HDFStore(snakemake.input[0], mode="r", locking=False) as store:
    for key in store.keys():
        meth_caller_dfs[key.strip("/")] = store[key]

# Heatmaps erstellen
chart = plot_heatmap_meth_callers(meth_caller_dfs)

# Chart speichern
chart.save(
    snakemake.output[0],
    embed_options={"actions": False},
    inline=False,
)
