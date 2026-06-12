import sys

import polars as pl

sys.stderr = open(snakemake.log[0], "w")
pl.Config.set_tbl_cols(20)


def parse_vcf(candidates):
    """Parses the VCF file and extracts relevant positions."""
    positions = []
    with open(candidates, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            chrom, pos = parts[0], int(parts[1])
            positions.append((chrom, pos))
    df = pl.DataFrame(positions, schema=["chrom", "pos"])
    return df


def parse_cov(file_path, strand, df):
    """Parses the coverage file and extracts coverage information."""
    cov_df = (
        pl.read_csv(
            file_path,
            separator="\t",
            has_header=False,
            schema={
                "chrom": pl.Utf8,
                "pos": pl.Int64,
                "end": pl.Int64,
                "coverage": pl.Float64,
            },
        )
        .select(["chrom", "pos", "coverage"])
        .with_columns(pl.col("coverage").alias(f"coverage_{strand}"))
        .drop("coverage")
    ).with_columns(
        pl.col("pos") + 1  # Convert to 1-based position
    )
    return df.join(
        cov_df,
        on=["chrom", "pos"],
        how="left",
    )


def ascii_to_methylation_expr(col: pl.Expr) -> pl.Expr:
    """
    Get the whole column of ascii values as input.
    Computes true methylation level according to https://github.com/seqan/seqan/blob/main/apps/mason2/README.mason_methylation
    """
    ascii_val = col.map_elements(ord, return_dtype=pl.Int32)
    # Skip '>' in order to have this character as indicator for chromosome
    threshold = ord(">")  # 62
    return (
        pl.when(ascii_val < threshold)
        .then((ascii_val - 33).cast(pl.Float64) / 80.0 * 100.0)
        .otherwise((ascii_val - 34).cast(pl.Float64) / 80.0 * 100.0)
    )


def parse_fasta(meth_file, df):
    """Parses a FASTA file and extracts methylation levels"""
    needed_positions = set(zip(df["chrom"].to_list(), df["pos"].to_list()))
    records_top = []
    records_bot = []
    current_chrom = None
    current_strand = None
    offset = 0

    with open(meth_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current_chrom, current_strand = line[1:].split("/")
                offset = 0
                continue

            if current_strand == "TOP":
                # Go through each character in the line and extract methylation levels
                for i, ascii_char in enumerate(line):
                    pos_1based = offset + i + 1
                    if (current_chrom, pos_1based) in needed_positions:
                        # Add to records if position is in candidates
                        records_top.append((current_chrom, pos_1based, ascii_char))

            elif current_strand == "BOT":
                for i, ascii_char in enumerate(line):
                    pos_0based = offset + i
                    if (current_chrom, pos_0based) in needed_positions:
                        records_bot.append((current_chrom, pos_0based, ascii_char))

            offset += len(line)

    lookup_top = pl.DataFrame(
        records_top,
        schema=["chrom", "pos", "ascii_top"],
        schema_overrides={"pos": pl.Int64},
    )
    lookup_bot = pl.DataFrame(
        records_bot,
        schema=["chrom", "pos", "ascii_bot"],
        schema_overrides={"pos": pl.Int64},
    )

    df = df.join(lookup_top, on=["chrom", "pos"], how="inner")
    df = df.join(lookup_bot, on=["chrom", "pos"], how="inner")
    df = df.with_columns(
        [
            ascii_to_methylation_expr(pl.col("ascii_top")).alias("meth_top"),
            ascii_to_methylation_expr(pl.col("ascii_bot")).alias("meth_bot"),
        ]
    ).drop(["ascii_top", "ascii_bot"])
    return df


def generate_bed(df):
    df = df.with_columns(
        (pl.col("coverage_TOP") + pl.col("coverage_BOT")).alias("total_coverage")
    ).with_columns(
        pl.when(pl.col("total_coverage") == 0)
        .then(0.0)
        .otherwise(
            (
                pl.col("meth_bot") * pl.col("coverage_BOT")
                + pl.col("meth_top") * pl.col("coverage_TOP")
            )
            / pl.col("total_coverage")
        )
        .alias("methylation_level")
    )
    return df


candidates = snakemake.input["candidates"]
meth_file = snakemake.input["methylation"]
cov_forward_file = snakemake.input["cov_forward"]
cov_reverse_file = snakemake.input["cov_reverse"]
output_file = snakemake.output[0]

df = parse_vcf(candidates)
df = parse_cov(cov_forward_file, "TOP", df)
df = parse_cov(cov_reverse_file, "BOT", df)
df = parse_fasta(meth_file, df)
df = generate_bed(df)
df.write_csv(output_file)
