import altair as alt
import pandas as pd
import polars as pl

# sys.stderr = open(snakemake.log[0], "w")

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)
pd.set_option("display.expand_frame_repr", False)
pd.set_option("display.max_rows", None)
pl.Config.set_tbl_cols(20)
pl.Config.set_tbl_rows(20)

sys.stderr = open(snakemake.log[0], "w")

alt.data_transformers.enable("vegafusion")

meth_caller_to_name = {
    "varlo_0.1": "Varlociraptor α = 0.1",
    "varlo_0.05": "Varlociraptor α = 0.05",
    "varlo_0.01": "Varlociraptor α = 0.01",
    "varlociraptor": "Varlociraptor",
    "bismark": "Bismark",
    "bsMap": "BSMAPz",
    "methylDackel": "MethylDackel",
    "modkit": "Modkit",
    "pb_CpG_tools": "Pb-CpG-tools",
    "bisSNP": "BisSNP",
}

tool_colors = {
    "Bismark": "#D81B60",
    "BSMAPz": "#1E88E5",
    "BisSNP": "#FFC107",
    "MethylDackel": "#f0700e",
    "Modkit": "#B42CEA",
    "Pb-CpG-tools": "#8D9279",
    "Varlociraptor α = 0.1": "#004D40",
    "Varlociraptor α = 0.05": "#126e5f",
    "Varlociraptor α = 0.01": "#05AA8F",
}


def bin_values(expr, bin_size):
    return ((expr / bin_size).round(0) * bin_size).cast(pl.Int64)


def plot_cov(df, caller, feature):

    all_bins_rep1 = list(range(0, df["rep1_cov_bin"].max() + 1))
    all_bins_rep2 = list(range(0, df["rep2_cov_bin"].max() + 1))
    ticks_axis_rep1 = list(range(0, df["rep1_cov_bin"].max() + 5, 5))
    ticks_axis_rep2 = list(range(0, df["rep2_cov_bin"].max() + 5, 5))
    x = (
        (df[feature] * df["count"]).sum() / df["count"].sum()
        if feature != "count"
        else df["count"].sum()
    )

    plot = (
        alt.Chart(df, title=f"{caller} - {feature}: {x:.2f}")
        .mark_rect()
        .encode(
            x=alt.X(
                "rep1_cov_bin:O",
                scale=alt.Scale(domain=all_bins_rep1),
                axis=alt.Axis(values=ticks_axis_rep1),
                title="Coverage replicate 1",
            ),
            y=alt.Y(
                "rep2_cov_bin:O",
                scale=alt.Scale(domain=list(reversed(all_bins_rep2))),
                axis=alt.Axis(values=ticks_axis_rep2),
                title="Coverage replicate 2",
            ),
            # color=alt.Color(
            #     f"{feature}:Q",
            #     scale=alt.Scale(
            #         domain=[-100, 0, 100],
            #         range=["red", "white", "blue"],
            #     ),
            # ),
            color=alt.Color(
                f"{feature}:Q",
            ),
            tooltip=[
                "rep1_cov_bin",
                "rep2_cov_bin",
                "mean_mae",
                "mean_mape",
                "count",
                "mean_methylation_rep1",
                "mean_methylation_rep2",
            ],
        )
        .properties(width=200, height=200)
    )
    return plot


def stratify_cov(df):
    # Calculate cutoffs ONCE for the entire dataset, not per caller
    rep1_cutoff = df.select(pl.col("rep1_cov_bin").quantile(quantile)).item()
    rep2_cutoff = df.select(pl.col("rep2_cov_bin").quantile(quantile)).item()

    # Now filter the entire df with these global cutoffs
    df_filtered = df.filter(
        (pl.col("rep1_cov_bin") <= rep1_cutoff)
        & (pl.col("rep2_cov_bin") <= rep2_cutoff)
    )

    mae_plots, mape_plots, count_plots = [], [], []
    all_dfs = {}
    for caller in meth_callers:
        # print(
        #     df_filtered.filter(pl.col("tool") == caller).filter(
        #         (pl.col("rep1_cov_bin") == 36) & (pl.col("rep2_cov_bin") == 57)
        #     )
        # )
        # Now aggregate per caller on the filtered data
        df_temp = (
            df_filtered.filter(pl.col("tool") == caller)
            .sort("rep1_cov_bin", "rep2_cov_bin")
            .group_by(["rep1_cov_bin", "rep2_cov_bin"])
            .agg(
                pl.col("mae").mean().alias("mean_mae"),
                pl.col("mape").mean().alias("mean_mape"),
                pl.col("methylation_rep1").mean().alias("mean_methylation_rep1"),
                pl.col("methylation_rep2").mean().alias("mean_methylation_rep2"),
                pl.len().alias("count"),
            )
            .to_pandas()
        )
        all_dfs[caller] = df_temp
        # Use the same ticks_axis for all callers

        mae_plots.append(plot_cov(df_temp, caller, "mean_mae"))
        mape_plots.append(plot_cov(df_temp, caller, "mean_mape"))
        count_plots.append(plot_cov(df_temp, caller, "count"))
    #######################################################
    # df_a = all_dfs["methylDackel"]
    # df_b = all_dfs["varlo_1.0"]
    # df_diff = df_a.merge(
    #     df_b, on=["rep1_cov_bin", "rep2_cov_bin"], suffixes=("_dackel", "_varlo")
    # )
    # df_diff["diff_mae"] = df_diff["mean_mae_dackel"] - df_diff["mean_mae_varlo"]
    # df_diff["diff_mape"] = df_diff["mean_mape_dackel"] - df_diff["mean_mape_varlo"]
    # df_diff["count"] = df_diff[["count_dackel", "count_varlo"]].min(axis=1)
    # print(df_diff)
    # diff_mae_plot = plot_cov(df_diff, "diff_mae", "diff_mae")
    # diff_mape_plot = plot_cov(df_diff, "diff_mape", "diff_mape")
    # chart = alt.vconcat(diff_mae_plot, diff_mape_plot).resolve_scale(
    #     color="independent"
    # )
    #######################################################

    mae_plot = alt.hconcat(*mae_plots)
    mape_plot = alt.hconcat(*mape_plots)
    count_plot = alt.hconcat(*count_plots)
    chart = alt.vconcat(mae_plot, mape_plot, count_plot).resolve_scale(
        color="independent"
    )
    return chart


def read_coverage(coverage_file: str, rep_name: str) -> pl.DataFrame:
    """Read and process coverage file, renaming coverage column appropriately."""
    return (
        pl.read_csv(
            coverage_file,
            separator="\t",
            has_header=False,
            new_columns=["chromosome", "pos_start", "pos_end", "coverage"],
        )
        .with_columns(
            (
                ((pl.col("pos_end") + pl.col("pos_start")) / 2)
                .cast(pl.Int64)
                .alias("position")
            )
        )
        .select(["chromosome", "position", "coverage"])
        .rename({"coverage": f"coverage_{rep_name}"})
    )


# Main execution
sample = snakemake.params["sample"]
plot_type = snakemake.params.get("plot_type")
coverage_bin_size = snakemake.params.get("bin_size", 1)
quantile = snakemake.params.get("quantile", 1.0)
meth_callers = snakemake.params["meth_callers"]


# Read and prepare methylation data
df = pl.read_parquet(snakemake.input["meth_data"])
df = df.drop([col for col in df.columns if "format" in col]).with_columns(
    pl.col("chromosome").cast(pl.Int64)
)
if sample != "all_samples":
    df = df.filter(pl.col("sample") == sample)
# Read and process coverage data
coverage_rep1 = read_coverage(snakemake.input["coverage_01"], "rep1")
coverage_rep2 = read_coverage(snakemake.input["coverage_02"], "rep2")
print(coverage_rep1)
df_cov_all = df.join(coverage_rep1, on=["chromosome", "position"], how="left").join(
    coverage_rep2, on=["chromosome", "position"], how="left"
)

df_long = (
    pl.concat(
        [
            df_cov_all.select(
                "chromosome",
                "position",
                pl.col(f"{caller}_methylation_rep1").alias("methylation_rep1"),
                pl.col(f"{caller}_methylation_rep2").alias("methylation_rep2"),
                "coverage_rep1",
                "coverage_rep2",
            ).with_columns(
                pl.lit(caller).alias("tool"),
            )
            for caller in meth_callers
        ]
    )
    .drop_nulls()
    .with_columns(
        bin_values(
            pl.min_horizontal(
                "methylation_rep1",
                "methylation_rep2",
            ),
            coverage_bin_size,
        ).alias("meth_bin"),
        pl.min_horizontal(
            "coverage_rep1",
            "coverage_rep2",
        ).alias("min_coverage"),
    )
)
df_min_coverage = (
    df_long.with_columns(
        bin_values(pl.col("min_coverage"), 5).alias("min_coverage_bin")
    )
    .group_by("tool", "min_coverage_bin")
    .agg(pl.len().alias("count"))
    .with_columns((pl.col("count") / pl.sum("count").over(["tool"])).alias("fraction"))
)
df_min_coverage = df_min_coverage.with_columns(
    pl.col("tool").replace(meth_caller_to_name).alias("tool_name")
)
color_domain = sorted(df_min_coverage["tool_name"].unique())
color_range = [tool_colors[t] for t in color_domain]
line_plot_min_cov_vs_count = (
    alt.Chart(df_min_coverage)
    .mark_line()
    .encode(
        x=alt.X("min_coverage_bin:Q", title="Min Coverage"),
        y=alt.Y("fraction:Q", title="Fraction"),
        color=alt.Color(
            "tool_name",
            title="Caller",
            scale=alt.Scale(
                domain=color_domain,
                range=color_range,
            ),
        ),
        strokeWidth=alt.value(1),
    )
)

# # Plot df_long as a line chart with color = tool, y-axis = mae and x-axis = meth_bin with altair
# df_meth = df_long.group_by("tool", "meth_bin").agg(
#     # pl.col("mae").mean(),
#     # pl.col("mape").mean(),
#     pl.len().alias("count"),
#     pl.col("min_coverage").mean().alias("mean_coverage"),
#     pl.col("min_coverage").median().alias("median_coverage"),
#     pl.col("min_coverage").quantile(0.01).alias("q01"),
#     pl.col("min_coverage").quantile(0.75).alias("q75"),
#     pl.col("min_coverage").min().alias("min_coverage"),
# )

# print(df_meth)

# # Plot histogram of min_coverage for meth_bin = 0, binned to size of 5, with one bar per tool, the bars per bin are next to each other
# df_meth_bin_0 = df_long.filter(pl.col("meth_bin") == 0)


# # Plot histogram of min_coverage for meth_bin = 5, binned to size of 5, with one bar per tool, the bars per bin are next to each other
# # df_meth_bin_5 = df_long.filter(pl.col("meth_bin") == 5)
# # Bin the min_coverage values using the existing bin_values function
# df_meth_bin = (
#     df_long.with_columns(
#         bin_values(pl.col("min_coverage"), 5).alias("min_coverage_bin")
#     )
#     .group_by("tool", "min_coverage_bin", "meth_bin")
#     .agg(pl.len().alias("count"))
#     .with_columns(
#         (pl.col("count") / pl.sum("count").over(["tool", "meth_bin"])).alias("fraction")
#     )
#     .sort("meth_bin")
# )
# df_meth_bin = df_meth_bin.filter(pl.col("meth_bin").is_in([0, 5, 100]))
# df_meth_bin = df_meth_bin.filter(pl.col("min_coverage_bin") < 100)

# histogram_meth_bin = (
#     alt.Chart(df_meth_bin)
#     .mark_line()
#     .encode(
#         x=alt.X("min_coverage_bin:O", title="Min Coverage (binned by 5)"),
#         xOffset="tool:N",
#         y="fraction:Q",
#         color="tool",
#         row="meth_bin:O",
#     )
# )


# min_coverage_by_meth_bin1 = (
#     alt.Chart(df_meth)
#     .mark_line()
#     .encode(
#         x="meth_bin",
#         y="q01",
#         color=alt.Color(
#             "tool:N",
#             title="Methylation caller",
#             scale=alt.Scale(range=colorblind_safe_palette.values()),
#             sort=meth_callers,
#         ),
#     )
# )

# min_coverage_by_meth_bin2 = (
#     alt.Chart(df_meth)
#     .mark_line()
#     .encode(
#         x="meth_bin",
#         y="q75",
#         color=alt.Color(
#             "tool:N",
#             title="Methylation caller",
#             scale=alt.Scale(range=colorblind_safe_palette.values()),
#             sort=meth_callers,
#         ),
#     )
# )
# min_coverage_by_meth_bin3 = (
#     alt.Chart(df_meth)
#     .mark_line()
#     .encode(
#         x="meth_bin",
#         y="median_coverage",
#         color=alt.Color(
#             "tool:N",
#             title="Methylation caller",
#             scale=alt.Scale(range=colorblind_safe_palette.values()),
#             sort=meth_callers,
#         ),
#     )
# )
# min_coverage_by_meth_bin4 = (
#     alt.Chart(df_meth)
#     .mark_line()
#     .encode(
#         x="meth_bin",
#         y="min_coverage",
#         color=alt.Color(
#             "tool:N",
#             title="Methylation caller",
#             scale=alt.Scale(range=colorblind_safe_palette.values()),
#             sort=meth_callers,
#         ),
#     )
# )

# min_coverage_by_meth_bin5 = (
#     alt.Chart(df_meth)
#     .mark_line()
#     .encode(
#         x="meth_bin",
#         y="mean_coverage",
#         color=alt.Color(
#             "tool:N",
#             title="Methylation caller",
#             scale=alt.Scale(range=colorblind_safe_palette.values()),
#             sort=meth_callers,
#         ),
#     )
# )

# line_chart = alt.concat(
#     # histogram_meth_bin_0,
#     line_plot_min_cov_vs_count,
#     histogram_meth_bin,
#     # min_coverage_by_meth_bin1,
#     # min_coverage_by_meth_bin2,
#     # min_coverage_by_meth_bin3,
#     # min_coverage_by_meth_bin4,
#     min_coverage_by_meth_bin5,
#     columns=2,
# )


line_plot_min_cov_vs_count.save(snakemake.output[0], scale_factor=2)
df_min_coverage.write_parquet(snakemake.output[1])
