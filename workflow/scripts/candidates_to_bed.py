import pysam

sys.stderr = open(snakemake.log[0], "w", buffering=1)

bcf_file = snakemake.input[0]

bcf = pysam.VariantFile(bcf_file)

with open(snakemake.output[0], "w") as bed_file:
    for record in bcf:
        chrom = record.chrom
        start = record.pos - 1  # BED ist 0-based
        end = record.pos  # VCF/BCF ist 1-based

        bed_file.write(f"{chrom}\t{start}\t{end}\n")
