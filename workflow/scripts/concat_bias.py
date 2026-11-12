import sys

import altair as alt
import pandas as pd

# Logging
sys.stderr = open(snakemake.log[0], "w")

# Same colors as in plot_bias.py
BIAS_COLORS = {
    "Strand Bias": "#D81B60",
    "Alt Locus Bias": "#1E88E5",
    "Read Orientation Bias": "#FFC107",
    "Read Position Bias": "#004D40",
    "Soft‑clipped Bias": "#8E44AD",
    "Haplotype Error": "#FF7F0E",
    "Fragment End Bias": "#6D4C41",
}

CATEGORIES = [
    "Bias both reps",
    "Bias, AF = 0",
    "Bias, AF ]0.0, 0.2]",
    "Bias, AF ]0.2, 0.4]",
    "Bias, AF ]0.4, 0.6]",
    "Bias, AF ]0.6, 0.8]",
    "Bias, AF ]0.8, 1.0]",
]


def bias_plots(df_long: pd.DataFrame, panel_order: list, alpha: str):
    """Create bias plots (one row per locus), one panel per platform/sample."""
    present_types = set(df_long["bias_type_label"])
    # Bias category plot
    bias_chart = (
        alt.Chart(df_long)
        .mark_bar()
        .encode(
            x=alt.X(
                "category:N",
                axis=alt.Axis(labelAngle=-45),
                title=None,
                scale=alt.Scale(domain=CATEGORIES),
            ),
            y=alt.Y("count():Q", title="Number of loci"),
            color=alt.Color(
                "bias_type_label:N",
                # Only the bias types present, so the legend has no unused entries
                scale=alt.Scale(
                    domain=[t for t in BIAS_COLORS if t in present_types],
                    range=[c for t, c in BIAS_COLORS.items() if t in present_types],
                ),
                title="Bias Type",
            ),
            tooltip=["category", "count()", "bias_type_label"],
            facet=alt.Facet(
                "platform_label:N", sort=panel_order, title=None, columns=3
            ),
        )
        .properties(title=f"ARTIFACT FDR α = {alpha}")
    )

    return bias_chart


def empty_plot(alpha):
    """Create an empty plot with a message."""
    chart = (
        alt.Chart(pd.DataFrame({"msg": [f"No bias data for ARTIFACT FDR α = {alpha}"]}))
        .mark_text(size=20)
        .encode(text="msg:N")
    )
    return chart


alphas = snakemake.params.alphas

# Read all bias dataframes and label each with its panel title
inputs = [
    (snakemake.input["illumina"], "Illumina"),
    (snakemake.input["pacbio"], "PacBio"),
    (snakemake.input["nanopore"], "Nanopore"),
]
panel_order = [label for _, label in inputs]

dfs = []
for path, label in inputs:
    df_panel = pd.read_parquet(path)
    df_panel["platform_label"] = label
    dfs.append(df_panel)
df = pd.concat(dfs, ignore_index=True)

# One row of platform panels per ARTIFACT FDR level α
charts = []
for alpha in alphas:
    df_alpha = df[df["fdr"] == alpha]
    if df_alpha.empty:
        charts.append(empty_plot(alpha))
    else:
        charts.append(bias_plots(df_alpha, panel_order, alpha))
chart = charts[0] if len(charts) == 1 else alt.vconcat(*charts)
chart.save(snakemake.output[0])
