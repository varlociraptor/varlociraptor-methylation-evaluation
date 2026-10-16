import sys

import altair as alt
import polars as pl

sys.stderr = open(snakemake.log[0], "w")

alt.data_transformers.enable("vegafusion")

# Same colors as in plot_discordance_vs_loci.py
TOOL_COLORS = {
    "Bismark": "#D81B60",
    "BSMAPz": "#1E88E5",
    "BisSNP": "#FFC107",
    "MethylDackel": "#f0700e",
    "Modkit": "#B42CEA",
    "Pb-CpG-tools": "#8D9279",
    "Varlociraptor α = 1.0": "#00362c",
    "Varlociraptor α = 0.5": "#005345",
    "Varlociraptor α = 0.1": "#00715f",
    "Varlociraptor α = 0.05": "#028d77",
    "Varlociraptor α = 0.01": "#05AA8F",
}

# Same names as in plot_discordance_vs_loci.py
METRIC_TO_NAME = {
    "mape": "Dᵣ",
    "mae": "Dₐ",
    # "binary_discordance": " Dₛ",
}


def plot_discordance_vs_loci(df: pl.DataFrame, sample_order: list[str]):
    color_domain = [
        tool for tool in TOOL_COLORS if tool in df["tool_name"].unique().to_list()
    ]
    color_range = [TOOL_COLORS[tool] for tool in color_domain]

    return (
        alt.Chart(df)
        .mark_line(point=alt.OverlayMarkDef(size=12), strokeWidth=1)
        .encode(
            x=alt.X("N:Q", title="Retained CpG loci (N)"),
            y=alt.Y("value:Q", title="Discordance", scale=alt.Scale(zero=False)),
            order=alt.Order("min_coverage:Q"),
            color=alt.Color(
                "tool_name:N",
                title="Caller",
                scale=alt.Scale(domain=color_domain, range=color_range),
            ),
            tooltip=[
                alt.Tooltip("tool_name:N", title="Caller"),
                alt.Tooltip("min_coverage:Q", title="Min coverage"),
                alt.Tooltip("N:Q", format=","),
                alt.Tooltip("value:Q", title="Discordance", format=".4f"),
            ],
        )
        .properties(width=200, height=200)
        .facet(
            row=alt.Row(
                "metric:N", title=None, sort=list(METRIC_TO_NAME.values())
            ),
            column=alt.Column("sample:N", title=None, sort=sample_order),
        )
        .resolve_scale(y="independent")
        .configure_header(labelFontSize=14, labelFontWeight="bold")
    )


platform_dfs = [
    pl.read_parquet(snakemake.input[platform])
    for platform in ["illumina", "pacbio", "nanopore"]
]
# Platforms side by side, samples within a platform sorted as in the single plots
sample_order = [
    sample for df in platform_dfs for sample in df["sample"].unique().sort()
]

chart = plot_discordance_vs_loci(pl.concat(platform_dfs), sample_order)
chart.save(snakemake.output[0], scale_factor=2)
