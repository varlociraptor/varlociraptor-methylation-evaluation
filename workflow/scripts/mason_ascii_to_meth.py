import re
import sys

import polars as pl

# Redirect standard error to snakemake log file
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


def parse_fasta(meth_file, df):
    """Parses a FASTA file and extracts methylation levels"""
    chrom_data = {}
    current_chrom = None
    current_strand = None

    with open(meth_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]
                current_chrom, current_strand = header.split("/")
                if current_chrom not in chrom_data:
                    chrom_data[current_chrom] = {"TOP": "", "BOT": ""}
            else:
                if current_chrom and current_strand:
                    chrom_data[current_chrom][current_strand] += line

    # Extract characters using the dictionary context directly
    df = df.with_columns(
        pl.struct(["chrom", "pos"])
        .map_elements(
            lambda x: (
                chrom_data[x["chrom"]]["TOP"][x["pos"] - 1]
                if x["chrom"] in chrom_data
                else None
            ),
            return_dtype=pl.Utf8,
        )
        .alias("ascii_top")
    ).with_columns(
        pl.struct(["chrom", "pos"])
        .map_elements(
            lambda x: (
                chrom_data[x["chrom"]]["BOT"][x["pos"]]
                if x["chrom"] in chrom_data
                else None
            ),
            return_dtype=pl.Utf8,
        )
        .alias("ascii_bot")
    )

    df = df.with_columns(
        [
            pl.col("ascii_top")
            .map_elements(ascii_to_methylation, return_dtype=pl.Float64)
            .alias("meth_top"),
            pl.col("ascii_bot")
            .map_elements(ascii_to_methylation, return_dtype=pl.Float64)
            .alias("meth_bot"),
        ]
    )

    return df.drop(["ascii_top", "ascii_bot"])


def ascii_to_methylation(char):
    """Compute true methylation level according to https://github.com/seqan/seqan/blob/main/apps/mason2/README.mason_methylation"""
    ascii_val = ord(char)
    if ascii_val < ord(">"):
        meth_level = ((ascii_val - ord("!")) / 80) * 100
    else:
        meth_level = ((ascii_val - ord("!") - 1) / 80) * 100
    return meth_level


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
    # """Generates a BED file with the collected methylation information."""
    # with open(output_file, "w") as out:
    #     for chrom, pos in candidate_positions:
    #         if chrom in meth_data and pos - 1 < len(meth_data[chrom]["BOT"]):
    #             coverage_top = coverages[(chrom, pos)]["TOP"]
    #             coverage_bot = coverages[(chrom, pos)]["BOT"]
    #             coverage = coverage_top + coverage_bot
    #             ascii_char_bot = meth_data[chrom]["BOT"][pos]
    #             ascii_char_top = meth_data[chrom]["TOP"][pos - 1]  # 1-based to 0-based
    #             meth_level = (
    #                 0
    #                 if coverage == 0
    #                 else (
    #                     ascii_to_methylation(ascii_char_bot) * coverage_bot
    #                     + ascii_to_methylation(ascii_char_top) * coverage_top
    #                 )
    #                 / coverage
    #             )
    #             if pos < 20:
    #                 print(
    #                     chrom,
    #                     pos,
    #                     ascii_char_bot,
    #                     ascii_char_top,
    #                     ascii_to_methylation(ascii_char_bot),
    #                     ascii_to_methylation(ascii_char_top),
    #                     coverage_bot,
    #                     coverage_top,
    #                     meth_level,
    #                 )

    #             out.write(
    #                 f"{chrom}\t{pos-1}\t{pos+1}\t{meth_level:.2f}\t{coverage}\t{ascii_char_bot}\t{ascii_char_top}\n"
    #             )


candidates = snakemake.input["candidates"]
meth_file = snakemake.input["methylation"]
cov_forward_file = snakemake.input["cov_forward"]
cov_reverse_file = snakemake.input["cov_reverse"]
output_file = snakemake.output[0]
df = parse_vcf(candidates)
df = parse_cov(cov_forward_file, "TOP", df)
df = parse_cov(cov_reverse_file, "BOT", df)
df = meth_data = parse_fasta(meth_file, df)
df = generate_bed(df)
df.write_csv(output_file)
