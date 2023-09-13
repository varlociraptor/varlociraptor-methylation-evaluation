import pickle

# Erweitern Sie die Wildcards und erstellen Sie eine Liste aller passenden Dateipfade

bedGraph_files = [snakemake.input[i] for i in range(len(snakemake.input))]

# Ein Dictionary zum Speichern der summierten Methylierungswerte pro Position
bedGrap_entry = {}

i = 0
# Iterieren Sie durch alle Bedgraph-Dateien
for file_path in bedGraph_files:
    with open(file_path, "r") as file:
        for line in file:
            if i % 50000 == 0:
                print(file_path, bedGraph_files, "\n")
            i+= 1
            parts = line.strip().split("\t")
            chrom, start, end, methylation, meth_reads, unmeth_reads = parts
            
            # Konvertieren Sie die Werte in Zahlen
            methylation = float(methylation)
            meth_reads = int(meth_reads)
            unmeth_reads = int(unmeth_reads)
            
            # Berechnen Sie den Durchschnitt für diese Position
            key = (chrom, start, end)
            if key not in bedGrap_entry:
                bedGrap_entry[key] = [methylation, meth_reads, unmeth_reads]
            else:
                bedGrap_entry[key][0] += meth_reads
                bedGrap_entry[key][1] += unmeth_reads

with open("bedGrap_entry.pkl", "wb") as outfile:
    pickle.dump(bedGrap_entry, outfile)

# Erstellen Sie eine neue Bedgraph-Datei für den Durchschnitt

output_file = snakemake.output[0]  # Name der Ausgabedatei
with open(output_file, "w") as outfile:
    for key, values in bedGrap_entry.items():
        if i % 1000 == 0:
            print(key, values)
        chrom, start, end = key
        meth_reads, unmeth_reads = values[0], values[1]
        avg_methylation = meth_reads / (meth_reads + unmeth_reads) if meth_reads + unmeth_reads != 0 else 0
        outfile.write(f"{chrom}\t{start}\t{end}\t{avg_methylation:.4f}\t{meth_reads:.0f}\t{unmeth_reads:.0f}\n")
        i += 1
