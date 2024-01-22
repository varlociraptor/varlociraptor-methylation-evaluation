import sys

# Einlesen der Eingabe- und Ausgabedateipfade aus den Snakemake-Argumenten
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Öffnen und Bearbeiten der Datei
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.startswith('>'):
            chrom_number = line[1:].strip()
            modified_line = f'>chr{chrom_number}\n'
            outfile.write(modified_line)
        else:
            outfile.write(line)
