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
alt.data_transformers.enable("vegafusion")


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
    title = f"Varlociraptor" if meth_caller == "varlo" else meth_caller.capitalize()

    # Compute MAPE if raw data is provided
    heatmap = (
        alt.Chart(
            df,
            title=alt.Title(
                title,
                subtitle=f"N = {df['count'].sum()} | D = {mapes[meth_caller]}",
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


def plot_histogram_cdf(meth_callers, meth_caller_dfs, cdf_dfs):
    df_hist_long = pd.concat(meth_caller_dfs.values(), axis=0, ignore_index=True)
    df_hist_long = df_hist_long.groupby(["dist_bin", "meth_caller"], as_index=False)[
        "rel_count"
    ].sum()
    # Make varlo always the same color
    domain = ["varlo"] + [c for c in meth_callers if c != "varlo"]

    df_cdf_long = pd.concat(cdf_dfs.values(), axis=0, ignore_index=True)

    # Histogram layer
    hist_chart = (
        alt.Chart(df_hist_long)
        .mark_point(size=60, filled=True)
        .encode(
            x=alt.X("dist_bin:Q", title="Distance (%)"),
            y=alt.Y("rel_count:Q", axis=alt.Axis(format="%")),
            color=alt.Color(
                "meth_caller:N", scale=alt.Scale(scheme="category10", domain=domain)
            ),
            tooltip=["dist_bin", "meth_caller", "rel_count"],
        )
    )

    # CDF layer (right y-axis)
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

    # Layer them
    combined_chart = alt.layer(cdf_chart, hist_chart).properties(
        width=700, height=400, title="CDF + Histogram of distances"
    )

    return combined_chart


def compute_replicate_counts(df_dict, bin_size, relative=False):
    """Compute replicate counts for each methylation caller"""
    protocols = snakemake.params["protocol"]
    meth_callers = snakemake.params["meth_callers"]

    # Allow heatmaps to be created per protocol or across multiple protocols by building the mean over all protocols
    print(df_dict.keys(), file=sys.stderr)
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
    print(df.to_string(), file=sys.stderr)
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
        print(meth_caller, temp_df.head(), file=sys.stderr)
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
        ).round(2)
        agg_cdf = agg_cdf.groupby("dist", as_index=False).agg(count=("dist", "size"))

        agg_cdf["rel_count"] = agg_cdf["count"] / agg_cdf["count"].sum() * 100
        agg_cdf = agg_cdf.sort_values("dist").reset_index(drop=True)
        agg_cdf["cdf"] = agg_cdf["count"].cumsum() / agg_cdf["count"].sum()
        agg_cdf["meth_caller"] = meth_caller
        print(agg_cdf.head(), file=sys.stderr)
        cdf_dfs[meth_caller] = agg_cdf[["dist", "cdf", "meth_caller"]]
    print(mapes, file=sys.stderr)
    
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


heatmap_plots = alt.hconcat(*heatmaps).resolve_scale(color="independent")

# distance_plots = plot_histogram_cdf(meth_callers, replicate_dfs, cdf_dfs)
# chart = alt.vconcat(heatmap_plots, distance_plots).resolve_scale(color="independent")

# Save final chart
heatmap_plots.save(
    snakemake.output[0],
    embed_options={"actions": False},
    inline=False,
)
