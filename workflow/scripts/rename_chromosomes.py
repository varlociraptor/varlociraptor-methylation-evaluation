import sys

def rename_chromosomes(input_file, output_file):
    print(input_file)
    print(output_file)
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        print(infile)
        for line in infile:
                print(line)
                if line.startswith('@'):  # Header line
                    outfile.write(line.replace('chr21', '21'))
                else:
                    fields = line.split('\t')
                    if fields[2] == 'chr21':
                        fields[2] = '21'
                    if fields[6] == 'chr21':
                        fields[6] = '21'
                    outfile.write('\t'.join(fields))

if __name__ == "__main__":
    input_file = snakemake.input.bam
    output_file = snakemake.output.renamed
    rename_chromosomes(input_file, output_file)
