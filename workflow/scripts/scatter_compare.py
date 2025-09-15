import pandas as pd
import altair as alt
import sys
import numpy as np

# Redirect error output to log file
sys.stderr = open(snakemake.log[0], "w")

# Show all columns during debugging
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)


# def plot_scatter_replicates(df_dict, corr_method):
#     charts_per_sample = []
#     sample_name = snakemake.params["protocol"]
#     print(df_dict.keys(), file=sys.stderr)
#     df = df_dict[sample_name]
#     sample_charts = []
#     for method in snakemake.params["methods"]:
#         # for method in ["varlo"]:
#         x_col = f"{method}_methylation_rep1"
#         y_col = f"{method}_methylation_rep2"

#         # if corr_method == "spearman":
#         #     df[f"{method}_rank_rep1"] = df[x_col].rank()
#         #     df[f"{method}_rank_rep2"] = df[y_col].rank()
#         #     x_col = f"{method}_rank_rep1"
#         #     y_col = f"{method}_rank_rep2"
#         df_temp = df.dropna(subset=[x_col, y_col])
#         # if x_col in df.columns and y_col in df.columns:
#         print(df_temp.head())
#         scatter = (
#             alt.Chart(df_temp)
#             .mark_rect()
#             .encode(
#                 x=alt.X(x_col, title=f"{method} Rep1", bin=alt.Bin(maxbins=10)),
#                 y=alt.Y(y_col, title=f"{method} Rep2", bin=alt.Bin(maxbins=10)),
#                 color=alt.Color(
#                     "count()", scale=alt.Scale(type="log", scheme="viridis")
#                 ),
#             )
#             .properties(
#                 title=f"Datapoints {len(df_temp)}",
#                 width=150,
#                 height=150,
#             )
#         )
#         print(scatter)
#         sample_charts.append(scatter)

#     if sample_charts:
#         charts_per_sample.append(
#             alt.hconcat(*sample_charts).properties(
#                 title=f"{sample_name} {corr_method} Correlation",
#             )
#         )
#     print(charts_per_sample)
#     return alt.vconcat(*charts_per_sample)


def bin_methylation(series, bin_size=10):
    return (series // bin_size) * bin_size


def plot_scatter_replicates(df_dict, corr_method):
    charts_per_sample = []
    sample_name = snakemake.params["protocol"]
    df = df_dict[sample_name]
    sample_charts = []
    for method in snakemake.params["methods"]:
        print(method)
        x_col = f"{method}_methylation_rep1"
        y_col = f"{method}_methylation_rep2"

        df_temp = df.dropna(subset=[x_col, y_col])

        df["rep1_bin"] = bin_methylation(df[x_col])
        df["rep2_bin"] = bin_methylation(df[y_col])

        agg = df.groupby(["rep1_bin", "rep2_bin"]).size().reset_index(name="count")

        scatter = (
            alt.Chart(agg)
            .mark_rect()
            .encode(
                x=alt.X(
                    "rep1_bin:O",
                    title=f"{method} Rep1",
                    sort=sorted(agg["rep1_bin"].unique(), reverse=False),  # absteigend
                ),
                y=alt.Y(
                    "rep2_bin:O",
                    title=f"{method} Rep2",
                    sort=sorted(agg["rep2_bin"].unique(), reverse=True),  # absteigend
                ),
                color=alt.Color(
                    "count:Q", scale=alt.Scale(type="log", scheme="viridis")
                ),
                tooltip=["rep1_bin", "rep2_bin", "count"],
            )
            .properties(
                title=f"Datapoints {len(df_temp)}",
                width=150,
                height=150,
            )
        )

        sample_charts.append(scatter)
    if sample_charts:
        charts_per_sample.append(
            alt.hconcat(*sample_charts).properties(
                title=f"{sample_name} {corr_method} Correlation",
            )
        )

    return alt.vconcat(*charts_per_sample)


replicate_dfs = {}
with pd.HDFStore(snakemake.input[0], mode="r", locking=False) as store:

    # with pd.HDFStore(snakemake.input[0]) as store:
    for key in store.keys():
        # key has leading '/', so remove it
        replicate_dfs[key.strip("/")] = store[key]


chart = plot_scatter_replicates(replicate_dfs, "pearson")

chart.save(
    snakemake.output[0],
    embed_options={"actions": False},
    inline=False,  # <- wichtig
)
