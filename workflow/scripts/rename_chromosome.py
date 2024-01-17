import sys

# Einlesen der Eingabe- und Ausgabedateipfade aus den Snakemake-Argumenten
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Öffnen und Bearbeiten der Datei
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        # Entfernen von 'chr' aus jeder Zeile
        modified_line = line.replace('chr', '')
        outfile.write(modified_line)