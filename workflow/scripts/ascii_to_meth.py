import re


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


def parse_fasta(meth_file):
    """Parses the ASCII FASTA file and extracts methylation levels."""
    methylation_data = {}
    with open(meth_file, "r") as f:
        current_seq = None
        for line in f:
            if line.startswith(">"):
                current_seq = line.strip().split("/")[0][1:]  # Extract chromosome ID
                methylation_data[current_seq] = ""
            else:
                methylation_data[current_seq] += line.strip()
    return methylation_data


def ascii_to_methylation(char):
    """Converts ASCII character to methylation level."""
    return (ord(char) - 33) * 1.25


def generate_bed(candidate_positions, meth_data, output_file):
    """Generates a BED file with methylation information."""
    with open(output_file, "w") as out:
        for chrom, pos in candidate_positions:
            if chrom in meth_data and pos - 1 < len(meth_data[chrom]):
                ascii_char = meth_data[chrom][pos]  # 1-based to 0-based
                meth_level = ascii_to_methylation(ascii_char)
                out.write(
                    f"{chrom}\t{pos-1}\t{pos+1}\t{meth_level:.2f}\t{ascii_char}\n"
                )


def main():
    candidates = snakemake.input["candidates"]
    meth_file = snakemake.input["methylation"]
    output_file = snakemake.output[0]

    candidate_positions = parse_vcf(candidates)
    meth_data = parse_fasta(meth_file)
    generate_bed(candidate_positions, meth_data, output_file)
    print(f"BED file created: {output_file}")


if __name__ == "__main__":
    main()
