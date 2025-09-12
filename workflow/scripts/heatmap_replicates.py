import pandas as pd
import altair as alt
import sys
import numpy as np
from functools import reduce

# Redirect error output to log file
sys.stderr = open(snakemake.log[0], "w")

# Show all columns during debugging
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)


def bin_methylation(series, bin_size):
    """Rundet Methylierungswerte auf die nächsten bin_size-Schritte und castet zu int"""
    return (np.round(series / bin_size) * bin_size).astype(int)


def plot_diff_heatmap(ref_tool, df_dict, bin_size):
    df1 = df_dict[ref_tool]
    df2 = df_dict["varlo"]

    # Merge counts
    df_diff = df1.merge(
        df2, on=["rep1_bin", "rep2_bin"], suffixes=(f"_{ref_tool}", "_varlo")
    )

    # Absolute Differenz
    df_diff["abs_diff"] = (df_diff[f"count_varlo"] - df_diff[f"count_{ref_tool}"]).abs()

    # log2 Fold Change
    df_diff["log2FC"] = np.where(
        (df_diff[f"count_varlo"] != 0) & (df_diff[f"count_{ref_tool}"] != 0),
        np.log2(df_diff[f"count_varlo"] / df_diff[f"count_{ref_tool}"]).abs(),
        0,
    )

    # Normalisierte Prozentwerte
    total_varlo = df_diff[f"count_varlo"].sum()
    total_ref_tool = df_diff[f"count_{ref_tool}"].sum()
    df_diff["rel_varlo"] = df_diff[f"count_varlo"] / total_varlo * 100
    df_diff["rel_ref_tool"] = df_diff[f"count_{ref_tool}"] / total_ref_tool * 100

    # Differenz in Prozent
    df_diff["rel_diff"] = df_diff["rel_varlo"] - df_diff["rel_ref_tool"]

    # Heatmap 1: Abs Difference
    diff_heatmap = (
        alt.Chart(
            df_diff, title=alt.Title(f"Abs Difference", subtitle=f"{ref_tool} vs varlo")
        )
        .mark_rect()
        .encode(
            x=alt.X("rep1_bin:O", sort=list(range(0, 101, bin_size))),
            y=alt.Y("rep2_bin:O", sort=list(range(100, -1, -bin_size))),
            color=alt.Color(
                "abs_diff:Q",
                scale=alt.Scale(
                    type="log", scheme="viridis", domain=[1, df_diff["abs_diff"].max()]
                ),
            ),
            tooltip=["rep1_bin", "rep2_bin", "abs_diff"],
        )
        .properties(width=200, height=200)
    )

    # Heatmap 2: Log2 Fold Change
    ratio_heatmap = (
        alt.Chart(
            df_diff,
            title=alt.Title(f"Log2 Fold Change", subtitle=f"{ref_tool} vs varlo"),
        )
        .mark_rect()
        .encode(
            x=alt.X("rep1_bin:O", sort=list(range(0, 101, bin_size))),
            y=alt.Y("rep2_bin:O", sort=list(range(100, -1, -bin_size))),
            color=alt.Color(
                "log2FC:Q",
                scale=alt.Scale(scheme="viridis", domain=[0, df_diff["log2FC"].max()]),
            ),
            tooltip=["rep1_bin", "rep2_bin", "log2FC"],
        )
        .properties(width=200, height=200)
    )

    # Heatmap 3: Relative Difference (Prozent)
    rel_heatmap = (
        alt.Chart(
            df_diff,
            title=alt.Title(f"Relative Difference", subtitle=f"{ref_tool} vs varlo"),
        )
        .mark_rect()
        .encode(
            x=alt.X("rep1_bin:O", sort=list(range(0, 101, bin_size))),
            y=alt.Y("rep2_bin:O", sort=list(range(100, -1, -bin_size))),
            color=alt.Color(
                "rel_diff:Q",
                scale=alt.Scale(
                    type="symlog",
                    scheme="redblue",  # divergierendes Schema
                    domain=[
                        -max(abs(df_diff["rel_diff"])),
                        max(abs(df_diff["rel_diff"])),
                    ],
                    domainMid=0,
                ),
            ),
            tooltip=["rep1_bin", "rep2_bin", "rel_diff"],
        )
        .properties(width=200, height=200)
    )

    return alt.vconcat(diff_heatmap, ratio_heatmap, rel_heatmap).resolve_scale(
        color="independent"
    )


def plot_heatmap_meth_callers(df_dict, bin_size):

    protocols = snakemake.params["protocol"]
    # We use this method to plot heatmaps for single protocols or common heatmaps over multiple protocols.
    # This depends on the snakemake parameter "protocol", which can be a string or a list of strings.
    if isinstance(protocols, str):
        df = df_dict[protocols]
    else:
        df_list = []
        for protocol in protocols:
            df_list.append(df_dict[protocol])
        long_df = pd.concat(df_list, axis=0, ignore_index=True)

        methyl_cols = [c for c in long_df.columns if "methylation" in c]
        df = long_df.groupby(["chromosome", "position"], as_index=False)[
            methyl_cols
        ].mean()
    meth_caller_dfs = {}
    method_plots = []
    distances_histo = pd.DataFrame()
    # Counts für jede meth_caller berechnen
    for meth_caller in snakemake.params["meth_callers"]:
        x_col = f"{meth_caller}_methylation_rep1"
        y_col = f"{meth_caller}_methylation_rep2"

        df_temp = df.dropna(
            subset=[x_col, y_col]
        ).copy()  # .copy(), um SettingWithCopyWarning zu vermeiden
        df_temp["rep1_bin"] = bin_methylation(df_temp[x_col], bin_size)
        df_temp["rep2_bin"] = bin_methylation(df_temp[y_col], bin_size)

        agg = (
            pd.crosstab(df_temp["rep1_bin"], df_temp["rep2_bin"])
            .stack()
            .reset_index(name="count")
        )
        ########### Histogram code

        def compute_dist(row):
            a, b = row["rep1_bin"], row["rep2_bin"]
            if max(a, b) == 0:
                return 0
            return abs((a - b) / max(a, b)) * 100

        # Is this measure the right way? We have a lot of distance 1.0, if one of the values is 0.
        agg["dist"] = agg.apply(compute_dist, axis=1)
        agg["dist_bin"] = (agg["dist"] / bin_size).round() * bin_size
        agg["dist_bin"] = agg["dist_bin"].astype(int)

        # Gruppieren nach Distanz-Bins
        dist_df = (
            agg.groupby("dist_bin", as_index=False)["count"]
            .sum()
            .assign(pct=lambda d: d["count"] / d["count"].sum() * 100)
            .sort_values("dist_bin")
            .rename(columns={"pct": meth_caller})
            .drop(columns="count")
        )

        if distances_histo.empty:
            distances_histo = dist_df
        else:
            distances_histo = distances_histo.merge(
                dist_df, on="dist_bin", how="outer"
            ).fillna(0)
        print(dist_df, file=sys.stderr)
        ############ Histogram code Ende

        meth_caller_dfs[meth_caller] = agg

        max_count = agg["count"].max()
        min_count = 1  # unterer Wert für log-Skala

        # Legend-Ticks setzen, inkl. Maximalwert
        ticks = [min_count, max_count]

        heatmap = (
            alt.Chart(
                agg,
                title=alt.Title(
                    f"{meth_caller} Heatmap", subtitle=f"Datapoints: {len(df_temp)}"
                ),
            )
            .mark_rect()
            .encode(
                x=alt.X(
                    "rep1_bin:O",
                    sort=list(range(0, 101, bin_size)),
                    title=f"{meth_caller} Rep1",
                ),
                y=alt.Y(
                    "rep2_bin:O",
                    sort=list(range(100, -1, -bin_size)),
                    title=f"{meth_caller} Rep2",
                ),
                color=alt.Color(
                    "count:Q",
                    scale=alt.Scale(
                        type="log", scheme="viridis", domain=[1, agg["count"].max()]
                    ),
                    legend=alt.Legend(
                        title="Count",
                        orient="right",
                        values=ticks,  # Ticks inkl. Maximalwert
                        format=",",
                    ),
                ),
                tooltip=["rep1_bin", "rep2_bin", "count"],
            )
            .properties(
                width=200,
                height=200,
            )
        )

        method_plots.append(heatmap)
    ## Histogram plotten
    print(distances_histo, file=sys.stderr)
    # DataFrame ins long format bringen
    df_long = distances_histo.melt(
        id_vars="dist_bin", var_name="meth_caller", value_name="count"
    )

    # Altair Chart (Balken nebeneinander statt gestapelt)
    distances_histo = (
        alt.Chart(df_long)
        .mark_bar()
        .encode(
            x=alt.X("dist_bin:O", title="Distanz"),  # Diskrete Bins auf der x-Achse
            y=alt.Y("count:Q", title="Count"),
            color=alt.Color("meth_caller:N", title="Meth Caller"),
            xOffset="meth_caller:N",  # 🔑 sorgt für nebeneinanderliegende Balken
            tooltip=["dist_bin", "meth_caller", "count"],
        )
        .properties(
            width=600, height=300, title="Histogramm der Distanzen pro Meth-Caller"
        )
    )

    return (
        alt.hconcat(*method_plots, distances_histo).resolve_scale(color="independent"),
        meth_caller_dfs,
    )


# Daten aus HDF5 einlesen
bin_size = snakemake.params["bin_size"]
meth_caller_dfs = {}
with pd.HDFStore(snakemake.input[0], mode="r", locking=False) as store:
    for key in store.keys():
        meth_caller_dfs[key.strip("/")] = store[key]

# Heatmaps erstellen
chart, meth_caller_dfs = plot_heatmap_meth_callers(meth_caller_dfs, bin_size)
diff_charts = []
# Differenz-Heatmaps erstellen
for meth_caller in snakemake.params["meth_callers"]:
    if meth_caller == "varlo":
        continue  # Differenz-Heatmap von varlo zu varlo macht keinen Sinn
    diff_chart = plot_diff_heatmap(meth_caller, meth_caller_dfs, bin_size)
    diff_charts.append(diff_chart)

diff_chart_plots = alt.hconcat(*diff_charts).resolve_scale(color="independent")
chart = alt.vconcat(chart, diff_chart_plots)
# Chart speichern
chart.save(
    snakemake.output[0],
    embed_options={"actions": False},
    inline=False,
)
