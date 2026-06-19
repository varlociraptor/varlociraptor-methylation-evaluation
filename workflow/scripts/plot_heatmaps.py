import sys

import altair as alt
import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
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
) -> alt.Chart:
    """Log-scaled heatmap for replicate methylation counts."""
    max_count = df["count"].max()
    ticks = list(np.logspace(0, np.log10(max_count), num=5).round().astype(int))
    mape = distances.loc[distances["meth_caller"] == meth_caller, "mape"].iloc[0]
    mae = distances.loc[distances["meth_caller"] == meth_caller, "mae"].iloc[0]
    binary_concordance = distances.loc[
        distances["meth_caller"] == meth_caller, "binary_concordance"
    ].iloc[0]
    heatmap = (
        alt.Chart(
            df,
            title=alt.Title(
                meth_caller_name,
                subtitle=f" N = {df['count'].sum():.0f} Dᵣ = {distance_df.filter(pl.col('meth_caller') == meth_caller)['mape'].item():.2f}%, Dₐ = {distance_df.filter(pl.col('meth_caller') == meth_caller)['mae'].item():.2f}%, Bc = {distance_df.filter(pl.col('meth_caller') == meth_caller)['binary_concordance'].item():.2f}",
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
fdr_levels = snakemake.params.get("fdr_levels", [])
meth_callers = combined_counts_df["meth_caller"].unique().tolist()

# Filter meth_callers to not include varlo_a with a not in fdr_levels
meth_callers = [
    m
    for m in meth_callers
    if not (m.split("_")[0] == "varlo" and float(m.split("_")[1]) not in fdr_levels)
]

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
heatmaps = [
    plot_heatmap(
        combined_counts_df[combined_counts_df["meth_caller"] == m],
        m,
        bin_size,
        distances,
        meth_caller_to_name.get(m, m),
    )
    for m in meth_callers
]
heatmap_plots = alt.hconcat(*heatmaps).resolve_scale(
    x="independent", y="independent", color="independent"
)


heatmap_plots.save(snakemake.output[0], embed_options={"actions": False}, inline=False)
