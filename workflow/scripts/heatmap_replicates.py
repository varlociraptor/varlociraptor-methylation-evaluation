import pandas as pd
import altair as alt
import sys
import numpy as np
from functools import reduce

# Redirect error output to log file
sys.stderr = open(snakemake.log[0], "w")

# Show all columns and many rows during debugging
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)


def bin_methylation(series, bin_size):
    """Round methylation values to nearest bin_size step and cast to int"""
    return (np.round(series / bin_size) * bin_size).astype(int)


def plot_count_heatmap(df, meth_caller, bin_size, mapes):
    """
    Create a log-scaled heatmap for replicate methylation values.
    Optionally compute MAPE if `raw_df` is provided (original methylation values).
    """
    max_count = df["count"].max()
    min_count = 1  # lower bound for log scale
    ticks = [min_count, max_count]  # legend ticks including maximum

    # Compute MAPE if raw data is provided
    heatmap = (
        alt.Chart(
            df,
            title=alt.Title(
                f"Absolute count {meth_caller} (log scaled)",
                subtitle=f"Datapoints: {df['count'].sum()} | MAPE: {mapes[meth_caller]}",
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
                    type="log", scheme="viridis", domain=[1, df["count"].max()]
                ),
                legend=alt.Legend(
                    title="Count",
                    orient="right",
                    values=ticks,
                    format=",",
                ),
            ),
            tooltip=["rep1_bin", "rep2_bin", "count"],
        )
        .properties(width=200, height=200)
    )
    return heatmap


def plot_histogram_cdf(meth_callers, meth_caller_dfs, cdf_dfs, bin_size=5):
    # all_dist_dfs = []
    # for meth_caller in meth_callers:
    #     df = meth_caller_dfs[meth_caller].copy()
    df_hist_long = pd.concat(meth_caller_dfs.values(), axis=0, ignore_index=True)
    print(df_hist_long, file=sys.stderr)
    df_hist_long = df_hist_long.groupby(["dist_bin", "meth_caller"], as_index=False)[
        "rel_count"
    ].sum()
    print(df_hist_long, file=sys.stderr)

    #     # CDF
    # print(meth_caller_dfs, file=sys.stderr)
    # print(cdf_dfs, file=sys.stderr)
    # # Merge histogram data
    # dist_df = all_dist_dfs[0]
    # for other in all_dist_dfs[1:]:
    #     dist_df = pd.merge(dist_df, other, on="dist_bin", how="outer").fillna(0)

    # df_hist_long = dist_df.melt(
    #     id_vars="dist_bin", var_name="meth_caller", value_name="rel_count"
    # )

    df_cdf_long = pd.concat(cdf_dfs.values(), axis=0, ignore_index=True)

    # Histogram layer
    hist_chart = (
        alt.Chart(df_hist_long)
        .mark_point(size=60, filled=True)
        .encode(
            x=alt.X("dist_bin:Q", title="Distance (%)"),
            y=alt.Y(
                "rel_count:Q", title="Relative number sites", axis=alt.Axis(format="%")
            ),
            color=alt.Color("meth_caller:N", scale=alt.Scale(scheme="category10")),
            tooltip=["dist_bin", "meth_caller", "rel_count"],
        )
    )

    # CDF layer (right y-axis)
    cdf_chart = (
        alt.Chart(df_cdf_long)
        .mark_line()
        .encode(
            x=alt.X("dist:Q"),
            y=alt.Y("cdf:Q", title="Cumulative fraction", axis=alt.Axis(format="%")),
            color=alt.Color("meth_caller:N", scale=alt.Scale(scheme="category10")),
            tooltip=["dist", "meth_caller", "cdf"],
        )
    )

    # Layer them
    combined_chart = alt.layer(cdf_chart, hist_chart).properties(
        width=700, height=400, title="Distance histogram and CDF combined"
    )

    return combined_chart


# def plot_dist_histogram(meth_callers, meth_caller_dfs, bin_size):
#     """Plot distance histograms between replicates across different methylation callers"""
#     all_dist_dfs = []

#     for meth_caller in meth_callers:
#         df = meth_caller_dfs[meth_caller]

#         df["dist"] = df.apply(
#             lambda row: (
#                 0
#                 if max(row["rep1_bin"], row["rep2_bin"]) == 0
#                 else abs(row["rep1_bin"] - row["rep2_bin"])
#                 / max(row["rep1_bin"], row["rep2_bin"])
#                 * 100
#             ),
#             axis=1,
#         )
#         df["dist_bin"] = (df["dist"] / bin_size).round() * bin_size
#         df["dist_bin"] = df["dist_bin"].astype(int)

#         dist_df = (
#             df.groupby("dist_bin", as_index=False)["rel_count"]
#             .sum()
#             .assign(pct=lambda d: d["rel_count"])
#             .sort_values("dist_bin")
#             .rename(columns={"pct": meth_caller})
#             .drop(columns="rel_count")
#         )
#         all_dist_dfs.append(dist_df)

#     # Merge distance distributions for all callers
#     dist_df = all_dist_dfs[0]
#     for other in all_dist_dfs[1:]:
#         dist_df = pd.merge(dist_df, other, on="dist_bin", how="outer").fillna(0)

#     # Convert to long format for plotting
#     df_long = dist_df.melt(
#         id_vars="dist_bin", var_name="meth_caller", value_name="rel_count"
#     )

#     line_chart = (
#         alt.Chart(df_long)
#         .mark_point(size=60, filled=True)
#         .encode(
#             x=alt.X(
#                 "dist_bin:O", title="Distance (%)", sort=list(range(0, 101, bin_size))
#             ),
#             y=alt.Y("rel_count:Q", title="Number of sites (%)"),
#             color=alt.Color(
#                 "meth_caller:N",
#                 title="Meth caller",
#                 scale=alt.Scale(scheme="category10"),
#             ),
#             tooltip=["dist_bin", "meth_caller", "rel_count"],
#         )
#         .properties(width=600, height=300, title="Distance histogram as line chart")
#     )
#     return line_chart


# def plot_cdf_distances(meth_callers, meth_caller_dfs):
#     """Plot cumulative distribution function (CDF) of replicate distances without binning"""
#     cdf_dfs = []

#     for meth_caller in meth_callers:
#         df = meth_caller_dfs[meth_caller].copy()

#         # Compute percentage distance between replicates
#         df["dist"] = df.apply(
#             lambda row: (
#                 0
#                 if max(row["rep1_bin"], row["rep2_bin"]) == 0
#                 else abs(row["rep1_bin"] - row["rep2_bin"])
#                 / max(row["rep1_bin"], row["rep2_bin"])
#                 * 100
#             ),
#             axis=1,
#         )

#         # Sort by distance
#         df = df.sort_values("dist").reset_index(drop=True)

#         # Compute weighted cumulative fraction
#         df["cdf"] = df["count"].cumsum() / df["count"].sum()
#         df["meth_caller"] = meth_caller

#         cdf_dfs.append(df[["dist", "cdf", "meth_caller"]])

#     # Combine all callers
#     df_cdf_long = pd.concat(cdf_dfs, axis=0, ignore_index=True)

#     # Plot CDF
#     cdf_chart = (
#         alt.Chart(df_cdf_long)
#         .mark_line()
#         .encode(
#             x=alt.X("dist:Q", title="Distance (%)"),
#             y=alt.Y("cdf:Q", title="Cumulative fraction"),
#             color=alt.Color(
#                 "meth_caller:N",
#                 title="Meth caller",
#                 scale=alt.Scale(scheme="category10"),
#             ),
#             tooltip=["dist", "meth_caller", "cdf"],
#         )
#         .properties(width=600, height=300, title="Cumulative distribution of distances")
#     )

#     return cdf_chart


# def plot_diff_heatmap(ref_tool, df_dict, bin_size):
#     """Plot absolute difference and log2 fold change heatmaps between reference and 'varlo' caller"""
#     df1 = df_dict[ref_tool]
#     df2 = df_dict["varlo"]

#     # Merge counts between reference and varlo
#     df_diff = df1.merge(
#         df2, on=["rep1_bin", "rep2_bin"], suffixes=(f"_{ref_tool}", "_varlo")
#     )

#     # Absolute difference
#     df_diff["diff"] = df_diff["rel_count_varlo"] - df_diff[f"rel_count_{ref_tool}"]

#     # Log2 fold change
#     df_diff["log2FC"] = np.where(
#         (df_diff["rel_count_varlo"] != 0) & (df_diff[f"rel_count_{ref_tool}"] != 0),
#         np.log2(df_diff["rel_count_varlo"] / df_diff[f"rel_count_{ref_tool}"]),
#         0,
#     )

#     diff_heatmap = (
#         alt.Chart(
#             df_diff,
#             title=alt.Title("Difference (log scaled)", subtitle=f"{ref_tool} vs varlo"),
#         )
#         .mark_rect()
#         .encode(
#             x=alt.X("rep1_bin:O", sort=list(range(0, 101, bin_size))),
#             y=alt.Y("rep2_bin:O", sort=list(range(100, -1, -bin_size))),
#             color=alt.Color(
#                 "diff:Q",
#                 scale=alt.Scale(
#                     type="symlog", scheme="redblue", domainMid=0
#                 ),  # force 0 to be mid (white)
#             ),
#             tooltip=["rep1_bin", "rep2_bin", "diff"],
#         )
#         .properties(width=200, height=200)
#     )

#     log2fc_heatmap = (
#         alt.Chart(
#             df_diff,
#             title=alt.Title(f"Log2 Fold Change", subtitle=f"{ref_tool} vs varlo"),
#         )
#         .mark_rect()
#         .encode(
#             x=alt.X("rep1_bin:O", sort=list(range(0, 101, bin_size))),
#             y=alt.Y("rep2_bin:O", sort=list(range(100, -1, -bin_size))),
#             color=alt.Color(
#                 "log2FC:Q",
#                 scale=alt.Scale(scheme="redblue", domainMid=0),
#             ),
#             tooltip=["rep1_bin", "rep2_bin", "log2FC"],
#         )
#         .properties(width=200, height=200)
#     )

#     return alt.vconcat(diff_heatmap, log2fc_heatmap).resolve_scale(color="independent")


def compute_replicate_counts(df_dict, bin_size, relative=False):
    """Compute replicate counts for each methylation caller"""
    protocols = snakemake.params["protocol"]
    meth_callers = snakemake.params["meth_callers"]

    # Allow heatmaps to be created per protocol or across multiple protocols
    if isinstance(protocols, str):
        df = df_dict[protocols]
    else:
        df_list = [df_dict[p] for p in protocols]
        long_df = pd.concat(df_list, axis=0, ignore_index=True)
        methyl_cols = [c for c in long_df.columns if "methylation" in c]
        df = long_df.groupby(["chromosome", "position"], as_index=False)[
            methyl_cols
        ].mean()

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

        df_temp = df.dropna(subset=[x_col, y_col]).copy()
        df_temp["rep1_bin"] = bin_methylation(df_temp[x_col], bin_size)
        df_temp["rep2_bin"] = bin_methylation(df_temp[y_col], bin_size)

        # We need all dist bin combos for our heatmap so we need to crosstab. Grouping is not enough!

        # Compute df with binned distances
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

        # Compute df with unbinned distances
        agg_cdf = df_temp.copy()
        agg_cdf["dist"] = agg_cdf.apply(
            lambda row: (
                0
                if max(row[x_col], row[y_col]) == 0
                else abs(row[x_col] - row[y_col]) / max(row[x_col], row[y_col]) * 100
            ),
            axis=1,
        )
        agg_cdf = agg_cdf.groupby("dist", as_index=False).agg(count=("dist", "size"))

        agg_cdf["rel_count"] = agg_cdf["count"] / agg_cdf["count"].sum() * 100
        agg_cdf = agg_cdf.sort_values("dist").reset_index(drop=True)
        agg_cdf["cdf"] = agg_cdf["count"].cumsum() / agg_cdf["count"].sum()
        agg_cdf["meth_caller"] = meth_caller
        cdf_dfs[meth_caller] = agg_cdf[["dist", "cdf", "meth_caller"]]
    return meth_caller_dfs, cdf_dfs, mapes


# Main execution ------------------------------------------------------

bin_size = snakemake.params["bin_size"]
meth_callers = snakemake.params["meth_callers"]
relative_counts = False

# Load data from HDF5
meth_caller_dfs = {}
with pd.HDFStore(snakemake.input[0], mode="r", locking=False) as store:
    for key in store.keys():
        meth_caller_dfs[key.strip("/")] = store[key]

# Compute replicate counts
replicate_dfs, cdf_dfs, mapes = compute_replicate_counts(
    meth_caller_dfs, bin_size, relative_counts
)

diff_charts = []
heatmaps = []
histogram_plots = []
# Create heatmaps and difference plots
for meth_caller in meth_callers:
    heatmaps.append(
        plot_count_heatmap(replicate_dfs[meth_caller], meth_caller, bin_size, mapes)
    )

distance_plots = plot_histogram_cdf(meth_callers, replicate_dfs, cdf_dfs)

# Combine plots
heatmap_plots = alt.hconcat(*heatmaps).resolve_scale(color="independent")
# diff_chart_plots = alt.hconcat(*diff_charts).resolve_scale(color="independent")
chart = alt.vconcat(heatmap_plots, distance_plots).resolve_scale(color="independent")

# Save final chart
chart.save(
    snakemake.output[0],
    embed_options={"actions": False},
    inline=False,
)
