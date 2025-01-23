import pandas as pd

def vcf_to_bed(vcf_file, bed_file):
    """
    Converts a VCF file containing CpG sites to a BED file.

    Args:
        vcf_file: Path to the input VCF file.
        bed_file: Path to the output BED file.
    """
    df = pd.read_csv(vcf_file, sep='\t', comment='#', names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    df['END'] = df['POS'] + 1
    df[['CHROM', 'POS', 'END']].to_csv(bed_file, sep='\t', index=False, header=False)

# Example usage
vcf_to_bed(snakemake.input[0], snakemake.output[0]) 