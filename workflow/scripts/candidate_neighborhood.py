candidate = snakemake.params["candidate"]
genome = snakemake.input["genome"]

# Die x-te Position aus der FASTA-Datei auslesen
with open(genome, "r") as fasta_file:
    header = fasta_file.readline()  # Die erste Zeile (Header) lesen und überspringen
    # Die restlichen Zeilen einlesen
    sequence = "".join([line.strip() for line in fasta_file])
    position_x = sequence[candidate - 1]

print(sequence[0:50])
left = candidate - 150 if candidate - 150 > 0 else 0
right = candidate + 150 if candidate + 150 < len(sequence) else len(sequence)


print(sequence[left:candidate-1], " ",
      sequence[candidate-1], " ", sequence[candidate:right])
