import sys

def rename_chromosomes(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('@'):  # Header line
                outfile.write(line.replace('chr', ''))
            else:
                fields = line.split('\t')
                if fields[2].startswith('chr'):
                    fields[2] = fields[2].removeprefix('chr')
                if fields[6].startswith('chr'):
                    fields[6] = fields[6].removeprefix('chr')
                if fields[9] != "*":
                    outfile.write('\t'.join(fields))

if __name__ == "__main__":
    input_file = snakemake.input.sam
    output_file = snakemake.output.renamed
    rename_chromosomes(input_file, output_file)
