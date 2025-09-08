import pandas as pd
import altair as alt
import sys
import numpy as np

# Redirect error output to log file
sys.stderr = open(snakemake.log[0], "w")

# Show all columns during debugging
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)


def plot_scatter_replicates(df_dict, corr_method):
    charts_per_sample = []
    sample_name = snakemake.params["protocol"]
    print(df_dict.keys(), file=sys.stderr)
    df = df_dict[sample_name]
    sample_charts = []
    for method in snakemake.params["methods"]:
        x_col = f"{method}_methylation_rep1"
        y_col = f"{method}_methylation_rep2"

        # if corr_method == "spearman":
        #     df[f"{method}_rank_rep1"] = df[x_col].rank()
        #     df[f"{method}_rank_rep2"] = df[y_col].rank()
        #     x_col = f"{method}_rank_rep1"
        #     y_col = f"{method}_rank_rep2"
        df_temp = df.dropna(subset=[x_col, y_col])
        # if x_col in df.columns and y_col in df.columns:
        scatter = (
            alt.Chart(df_temp)
            .mark_circle(size=10, opacity=0.4)
            .encode(
                x=alt.X(x_col, title=f"{method} Rep1"),
                y=alt.Y(y_col, title=f"{method} Rep2"),
                color=alt.Color('count()',scale=alt.Scale(scheme='viridis')),
                # tooltip=["chromosome", "position", x_col, y_col],
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
# def plot_heatmap_replicates(df_dict, corr_method):
#     charts_per_sample = []
#     sample_name = snakemake.params["protocol"]
#     df = df_dict[sample_name]
#     sample_charts = []
#     # print("df_temp:", file=sys.stderr)

#     # print(df_dict.hea/d(), file=sys.stderr)
#     for method in snakemake.params["methods"]:
#         x_col = f"{method}_methylation_rep1"
#         y_col = f"{method}_methylation_rep2"

#         df_temp = df.dropna(subset=[x_col, y_col])
#         print("Length df_temp:", file=sys.stderr)
#         print(df_temp.head(), file=sys.stderr)

#         # 🟢 Voraggregieren in 100x100 bins mit numpy
#         H, xedges, yedges = np.histogram2d(
#             df_temp[x_col], df_temp[y_col],
#             bins=100, range=[[0, 100], [0, 100]]  # anpassen falls Wertebereich ≠ 0..1
#         )
#         print("Heatmap df_temp:", len(df_temp), file=sys.stderr)

#         heatmap_df = pd.DataFrame({
#             "x": np.repeat((xedges[:-1] + xedges[1:]) / 2, 100),
#             "y": np.tile((yedges[:-1] + yedges[1:]) / 2, 100),
#             "count": H.flatten()
#         })
#         print(heatmap_df, file=sys.stderr)


#         heatmap = (
#             alt.Chart(heatmap_df)
#             .mark_rect()
#             .encode(
#                 x=alt.X("x:Q", title=f"{method} Rep1").scale(type="log"),
#                 y=alt.Y("y:Q", title=f"{method} Rep2").scale(type="log"),
#                 color=alt.Color("count:Q", scale=alt.Scale(scheme="viridis"))
#             )
#             .properties(
#                 title=f"Datapoints {len(df_temp)}",
#                 width=150,
#                 height=150,
#             )
#         )

#         sample_charts.append(heatmap)

#     if sample_charts:
#         charts_per_sample.append(
#             alt.hconcat(*sample_charts).properties(
#                 title=f"{sample_name} {corr_method} Heatmap",
#             )
#         )

#     return alt.vconcat(*charts_per_sample)



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
