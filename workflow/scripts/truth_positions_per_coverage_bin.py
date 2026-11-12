import sys

import polars as pl

sys.stderr = open(snakemake.log[0], "w")

coverage_bins = snakemake.params.coverage_bins
i = int(snakemake.wildcards.bin)

truth_df = pl.read_csv(snakemake.input[0], schema_overrides={"chrom": pl.Utf8})

# Bin i covers [coverage_bins[i], coverage_bins[i + 1]), the last bin has no upper bound
in_bin = pl.col("total_coverage") >= coverage_bins[i]
if i + 1 < len(coverage_bins):
    in_bin &= pl.col("total_coverage") < coverage_bins[i + 1]

# Tab-separated CHROM and POS without header, as expected by bcftools view -T
truth_df.filter(in_bin).select("chrom", "pos").write_csv(
    snakemake.output[0], separator="\t", include_header=False
)
