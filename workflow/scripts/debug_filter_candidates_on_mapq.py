import pysam
import cyvcf2

# Redirect standard error to snakemake log file
sys.stderr = open(snakemake.log[0], "w")


def filter_variants(
    bcf_file, bam_file, output_bcf, mapq_threshold=60, min_fraction=0.5
):
    """
    Filter variants from a BCF file based on mapping quality (MAPQ) of aligned reads.
    """
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Open the BCF file
    vcf_reader = cyvcf2.VCF(bcf_file)
    vcf_writer = cyvcf2.Writer(output_bcf, vcf_reader)

    for record in vcf_reader:
        chrom = record.CHROM
        pos = record.POS

        # Fetch reads at the given position
        reads = bam.fetch(chrom, pos - 1, pos)

        mapq_values = [read.mapping_quality for read in reads if not read.is_unmapped]
        if not mapq_values:
            continue  # No reads found

        # Calculate the fraction of reads with MAPQ < threshold
        low_mapq_count = sum(1 for q in mapq_values if q < mapq_threshold)
        fraction_low_mapq = low_mapq_count / len(mapq_values)

        if fraction_low_mapq >= min_fraction:
            vcf_writer.write_record(record)

    # Close files
    bam.close()
    vcf_writer.close()


# Example call

filter_variants(
    snakemake.input["candidates"], snakemake.input["alignment"], snakemake.output[0]
)
