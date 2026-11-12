import sys

import altair as alt
import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)
alt.data_transformers.enable("vegafusion")

titles = {
    "np_methylSeq": "Nanopore and MethylSeq",
    "pb_methylSeq": "PacBio and MethylSeq",
    "np_pb": "Nanopore and PacBio",
}


def plot_heatmap(
    df: pd.DataFrame,
    meth_caller: str,
    bin_size: int,
    distances: dict,
    meth_caller_name: str,
    max_count: int,
) -> alt.Chart:
    """Log-scaled heatmap for replicate methylation counts."""
    ticks = list(np.logspace(0, np.log10(max_count), num=5).round().astype(int))
    mape = distances.loc[distances["meth_caller"] == meth_caller, "mape"].iloc[0]
    mae = distances.loc[distances["meth_caller"] == meth_caller, "mae"].iloc[0]
    binary_discordance = distances.loc[
        distances["meth_caller"] == meth_caller, "binary_discordance"
    ].iloc[0]
    number_points = df["count"].sum() // 1000
    # Multi-sample: comparison as title; in the report additionally the caller
    # (incl. FDR level) as second line, in the paper heatmaps only the comparison
    comparison = snakemake.output[0].split("/")[-3]
    if comparison in titles:
        title = (
            [titles[comparison], meth_caller_name]
            if snakemake.params["alpha_in_title"]
            else titles[comparison]
        )
    else:
        title = meth_caller_name

    heatmap = (
        alt.Chart(
            df,
            title=alt.TitleParams(
                text=title,
                subtitle=[
                    f"N = {number_points}k | Dᵣ = {mape:.2f} | Dₐ = {mae:.2f}",
                ],
                subtitleFontSize=9,
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


combined_counts_df = pd.read_parquet(snakemake.input["df"], engine="pyarrow")
distances = pd.read_parquet(snakemake.input["distances"], engine="pyarrow")
bin_size = snakemake.params["bin_size"]

available_callers = set(combined_counts_df["meth_caller"].unique())
meth_callers = [m for m in snakemake.params["callers"] if m in available_callers]
missing_callers = [m for m in snakemake.params["callers"] if m not in available_callers]


plot_type = snakemake.params.get("plot_type")
meth_caller_to_name = {
    "bismark": "Bismark",
    "bsMap": "BSMAPz",
    "bisSNP": "BisSNP",
    "methylDackel": "MethylDackel",
    "modkit": "Modkit",
    "pb_CpG_tools": "pb-CpG-tools",
}
for m in meth_callers:
    if m.startswith("varlo_"):
        alpha = m.split("_")[1]
        meth_caller_to_name[m] = f"Varlociraptor α = {alpha}"
max_count = combined_counts_df.loc[
    combined_counts_df["meth_caller"].isin(meth_callers), "count"
].max()
heatmaps = [
    plot_heatmap(
        combined_counts_df[combined_counts_df["meth_caller"] == m],
        m,
        bin_size,
        distances,
        meth_caller_to_name.get(m, m),
        max_count
    )
    for m in meth_callers
]
heatmap_plots = alt.concat(*heatmaps, columns=4).resolve_scale(
    x="independent", y="independent", color="shared"
).properties()


heatmap_plots.save(snakemake.output[0], embed_options={"actions": False}, inline=False)
