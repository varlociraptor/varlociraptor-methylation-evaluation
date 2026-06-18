from pathlib import Path

import polars as pl

sys.stderr = open(snakemake.log[0], "w")


dfs = []

for f in snakemake.input:
    df = pl.read_csv(
        f,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end", "coverage"],
    ).with_columns(pl.col("coverage").cast(pl.Float64))
    dfs.append(df)

merged = (
    pl.concat(dfs)
    .group_by(["chrom", "start", "end"])
    .agg(pl.col("coverage").mean().alias("coverage"))
    .sort(["chrom", "start"])
)

merged.write_csv(snakemake.output[0], separator="\t", include_header=False)
