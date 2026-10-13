import re
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
            # Skip calls without AFD, as in the PIT calibration
            format_dict = dict(zip(parts[8].split(":"), parts[9].split(":")))
            afd = format_dict.get("AFD")
            if afd is None or afd == ".":
                continue
            records.append((parts[0], int(parts[1]), phred_to_prob(prob_high)))
    return pl.DataFrame(
        records, schema=["chrom", "pos", "prob_high"]).with_columns((1 - pl.col("prob_high")).alias("prob_not_high"))


def coverage_label(bins, i):
    if i == len(bins) - 1:
        return f"≥ {bins[i]}"
    return f"{bins[i]}–{bins[i + 1] - 1}"


truth_df = pl.read_csv(snakemake.input.truth, schema_overrides={"chrom": pl.Utf8})

high_threshold = snakemake.params.high_threshold
coverage_bins = snakemake.params.coverage_bins
alphas = snakemake.params.fdr_levels


coverage_labels = [coverage_label(coverage_bins, i) for i in range(len(coverage_bins))]
truth_df = truth_df.with_columns(
    # Boolean if position is considered highly methylated in truth
    (pl.col("methylation_level") > high_threshold).alias("true_high"),
    # Assign a coverage bin to each single coverage value
    pl.col("total_coverage")
    .cut(coverage_bins[1:], labels=coverage_labels, left_closed=True)
    .cast(pl.Utf8)
    .alias("coverage"),
)
# Number of sites per coverage bin before FDR control
num_sites = {}
for unfiltered in snakemake.input.unfiltered:
    # File names look like bin_{bin}.vcf
    bin_index = int(re.search(r"bin_(\d+)\.vcf$", unfiltered).group(1))
    num_sites[bin_index] = len(parse_calls(unfiltered))

rows = []

for filtered in snakemake.input.filtered:
    # File names look like bin_{bin}_fdr_{alpha}.vcf
    bin_index, alpha = re.search(r"bin_(\d+)_fdr_([\d.]+)\.vcf$", filtered).groups()
    bin_index = int(bin_index)
    # Calls selected by varlociraptor's global FDR control within this coverage bin
    selected = parse_calls(filtered).join(truth_df, on=["chrom", "pos"], how="inner")
    num_calls = len(selected)
    num_false = num_calls - selected["true_high"].sum()

    rows.append(
        {
            "bin": bin_index,
            "coverage": coverage_labels[bin_index],
            "num_sites": num_sites[bin_index],
            "fdr_threshold": float(alpha),
            "num_calls": num_calls,
            "num_false_discoveries": num_false,
            # Mean posterior probability of not being HIGH over the selected calls
            "expected_fdr": selected["prob_not_high"].mean() if num_calls else None,
            "true_fdr": num_false / num_calls if num_calls else None,
        }
    )

# Order by coverage bin, as the legend takes the bin order from the rows
result_df = pl.DataFrame(rows).sort("bin", "fdr_threshold").drop("bin")
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
    plot_df = result_df.with_columns(
        pl.format("{} (n={})", "coverage", "num_sites").alias("coverage")
    )
    # Rows are ordered by coverage bin
    legend_labels = plot_df["coverage"].unique(maintain_order=True).to_list()
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
                scale=alt.Scale(domain=legend_labels, range=coverage_colors),
                title="Coverage",
            ),
        )
    )
    return alt.layer(diagonal, curves).properties(width=250, height=250)


chart = plot_fdr()
chart.save(snakemake.output.plot)
