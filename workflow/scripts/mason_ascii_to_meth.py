from collections import defaultdict
import re


# Redirect standard error to snakemake log file
sys.stderr = open(snakemake.log[0], "w")


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
    return positions


def parse_cov(file_path, strand, coverages):
    """Parses the coverage file and extracts coverage information."""
    with open(file_path, "r") as file:

        for line in file:
            parts = line.strip().split("\t")
            chrom = parts[0]
            position = int(parts[1])
            coverage = float(parts[3])

            coverages[(chrom, position)][strand] = coverage
    return coverages


def parse_fasta(meth_file):
    """Parses the ASCII FASTA file and extracts methylation levels."""
    with open(meth_file, "r") as f:
        fasta = f.read().strip()

        fasta_parts = re.split(r">|/TOP\n|/BOT\n", fasta)
        chrom = fasta_parts[1]
        methylation_data = defaultdict(lambda: defaultdict(str))
        methylation_data[chrom]["TOP"] = fasta_parts[2].replace("\n", "")
        methylation_data[chrom]["BOT"] = fasta_parts[4].replace("\n", "")
        return methylation_data


def ascii_to_methylation(char):
    """Compute true methylation level according to https://github.com/seqan/seqan/blob/main/apps/mason2/README.mason_methylation"""
    ascii_val = ord(char)
    if ascii_val < ord(">"):
        L = ((ascii_val - ord("!")) / 80) * 100
    else:
        L = ((ascii_val - ord("!") - 1) / 80) * 100
    return L


def generate_bed(candidate_positions, meth_data, coverages, output_file):
    """Generates a BED file with the collected methylation information."""
    with open(output_file, "w") as out:
        for chrom, pos in candidate_positions:
            if chrom in meth_data and pos - 1 < len(meth_data[chrom]["BOT"]):
                coverage_top = coverages[(chrom, pos)]["TOP"]
                coverage_bot = coverages[(chrom, pos)]["BOT"]
                coverage = coverage_top + coverage_bot
                ascii_char_bot = meth_data[chrom]["BOT"][pos]
                ascii_char_top = meth_data[chrom]["TOP"][pos - 1]  # 1-based to 0-based
                meth_level = (
                    0
                    if coverage == 0
                    else (
                        ascii_to_methylation(ascii_char_bot) * coverage_bot
                        + ascii_to_methylation(ascii_char_top) * coverage_top
                    )
                    / coverage
                )

                out.write(
                    f"{chrom}\t{pos-1}\t{pos+1}\t{meth_level:.2f}\t{coverage}\t{ascii_char_bot}\t{ascii_char_top}\n"
                )


candidates = snakemake.input["candidates"]
meth_file = snakemake.input["methylation"]
cov_forward_file = snakemake.input["cov_forward"]
cov_reverse_file = snakemake.input["cov_reverse"]
output_file = snakemake.output[0]
candidate_positions = parse_vcf(candidates)
cov_forward = parse_cov(
    cov_forward_file, "TOP", defaultdict(lambda: defaultdict(float))
)
coverages = parse_cov(cov_reverse_file, "BOT", cov_forward)
meth_data = parse_fasta(meth_file)
generate_bed(candidate_positions, meth_data, coverages, output_file)
