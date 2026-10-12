import sys

import altair as alt
import numpy as np
import polars as pl

sys.stderr = open(snakemake.log[0], "w")


def phred_to_prob(score: str) -> float:
    return 10 ** (-float(score) / 10)


def parse_calls(filepath: str) -> pl.DataFrame:
    """Extract the posterior probability of the HIGH event for every call."""
    records = []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[4] != "<METH>":
                continue
            info_dict = dict(
                item.split("=", 1) for item in parts[7].split(";") if "=" in item
            )
            prob_high = info_dict.get("PROB_HIGH")
            if prob_high is None or prob_high == ".":
                continue
            records.append((parts[0], int(parts[1]), phred_to_prob(prob_high)))
    return pl.DataFrame(
        records, schema=["chrom", "pos", "prob_high"]).with_columns((1 - pl.col("prob_high")).alias("prob_not_high"))


def coverage_label(bins, i):
    if i == len(bins) - 1:
        return f"≥ {bins[i]}"
    return f"{bins[i]}–{bins[i + 1] - 1}"


def select_local(df: pl.DataFrame, alpha: float) -> pl.DataFrame:
    """Keep every site whose posterior probability of not being HIGH is at most alpha."""

    return df.filter(pl.col("prob_high") >= 1 - alpha)



truth_df = pl.read_csv(snakemake.input.truth, schema_overrides={"chrom": pl.Utf8})
calls_df = parse_calls(snakemake.input.calls)

high_threshold = snakemake.params.high_threshold
coverage_bins = snakemake.params.coverage_bins
alphas = snakemake.params.fdr_levels


coverage_labels = [coverage_label(coverage_bins, i) for i in range(len(coverage_bins))]
df = calls_df.join(truth_df, on=["chrom", "pos"], how="inner").with_columns(
    # Boolean if position is considered highly methylated in truth
    (pl.col("methylation_level") > high_threshold).alias("true_high"),
    # Assign a coverage bin to each single coverage value
    pl.col("total_coverage")
    .cut(coverage_bins[1:], labels=coverage_labels, left_closed=True)
    .cast(pl.Utf8)
    .alias("coverage"),
)
# Divide df into smaller dfs for each coverage bin
strata = [("all", df)] + [
    (label, df.filter(pl.col("coverage") == label)) for label in coverage_labels
]

rows = []

for stratum, stratum_df in strata:
    for alpha in alphas:
        selected = select_local(stratum_df, alpha)
        num_calls = len(selected)
        num_false = num_calls - selected["true_high"].sum()

        rows.append(
            {
                "coverage": stratum,
                "num_sites": len(stratum_df),
                "fdr_threshold": float(alpha),
                "num_calls": num_calls,
                "num_false_discoveries": num_false,
                "true_fdr": num_false / num_calls if num_calls else None,
            }
        )

result_df = pl.DataFrame(rows)
result_df.write_csv(snakemake.output.table)

coverage_colors = ["#a8e6d9", "#7fd8c6", "#05AA8F", "#126e5f", "#012b24"]


fdr_levels = snakemake.params.fdr_levels
diagonal = (
    alt.Chart(pl.DataFrame({"x": [0, 1]}))
    .mark_line(strokeDash=[4, 4], color="#888888", strokeWidth=1)
    .encode(x=alt.X("x:Q"), y=alt.Y("x:Q"))
)


ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]


def plot_fdr():
    plot_df = result_df.filter(pl.col("coverage") != "all")
    curves = (
        alt.Chart(plot_df)
        .mark_line(point=alt.OverlayMarkDef(size=30), strokeWidth=2)
        .encode(
            x=alt.X(
                "fdr_threshold:Q",
                scale=alt.Scale(domain=[0, 1]),
                axis=alt.Axis(values=ticks),
                title="FDR threshold (α)",
            ),
            y=alt.Y(
                "true_fdr:Q",
                scale=alt.Scale(domain=[0, 1]),
                axis=alt.Axis(values=ticks),
                title="true FDR",
            ),
            color=alt.Color(
                "coverage:N",
                scale=alt.Scale(domain=coverage_labels, range=coverage_colors),
                title="Coverage",
            ),
        )
    )
    return alt.layer(diagonal, curves).properties(width=250, height=250)


chart = plot_fdr()
chart.save(snakemake.output.plot)
