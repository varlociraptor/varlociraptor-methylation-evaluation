import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


def vcf_to_bed(vcf_file, bed_file):
    """Converts a VCF file containing CpG sites to a BED file."""
    df = pd.read_csv(
        vcf_file,
        sep="\t",
        comment="#",
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
    )
    df["END"] = df["POS"] + 1
    df[["CHROM", "POS", "END"]].to_csv(bed_file, sep="\t", index=False, header=False)


vcf_to_bed(snakemake.input[0], snakemake.output[0])
