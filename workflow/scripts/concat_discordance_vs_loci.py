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
    # One Varlociraptor curve over all α
    "Varlociraptor": "#028d77",
}

# Same names as in plot_discordance_vs_loci.py
METRIC_TO_NAME = {
    "mape": "Relative discordance",
    "mae": "Absolute discordance",
    # "binary_discordance": " Dₛ",
}

# Points based on fewer loci are too noisy and are not shown
MIN_LOCI = 5000


def plot_discordance_vs_loci(df: pl.DataFrame, sample_order: list[str]):
    df = df.filter(pl.col("N") > MIN_LOCI).with_columns(
        # Rightmost point of each curve (largest N) is drawn bigger
        pl.when(pl.col("N") == pl.col("N").max().over("sample", "metric", "tool_name"))
        .then(40)
        .otherwise(20)
        .alias("point_size")
    )
    color_domain = [
        tool for tool in TOOL_COLORS if tool in df["tool_name"].unique().to_list()
    ]
    color_range = [TOOL_COLORS[tool] for tool in color_domain]
    base = alt.Chart().encode(
        x=alt.X("N:Q", title="Retained CpG loci (N)"),
        y=alt.Y("value:Q", title=None, scale=alt.Scale(zero=False)),
        color=alt.Color(
            "tool_name:N",
            title="Caller",
            scale=alt.Scale(domain=color_domain, range=color_range),
        ),
    )
    line = base.mark_line(strokeWidth=1).encode(order=alt.Order("N:Q"))
    points = base.mark_point(filled=True, opacity=1).encode(
        size=alt.Size("point_size:Q", scale=None, legend=None)
    )

    return (
        alt.layer(line, points, data=df)
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
# vl-convert rejects scale != 1 for PDF (vega_to_pdf path, taken when vegafusion pre-compiles the spec)
chart.save(
    snakemake.output[0],
    scale_factor=1 if snakemake.output[0].endswith(".pdf") else 2,
)
