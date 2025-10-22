import sys

# Rename chromosome names in a FASTA file by adding 'chr' prefix.

sys.stderr = open(snakemake.log[0], "w")

input_file = snakemake.input[0]
output_file = snakemake.output[0]

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            chrom_number = line[1:].strip()
            modified_line = f">chr{chrom_number}\n"
            outfile.write(modified_line)
        else:
            outfile.write(line)
