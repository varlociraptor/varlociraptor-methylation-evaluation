#  Redirect standard error to snakemake log file
sys.stderr = open(snakemake.log[0], "w")


# Gather input bedGraph files from Snakemake
bedGraph_files = [file for file in snakemake.input]

# Dictionary to store methylation data
bedGraph_entry = {}

# Process each input file
for file_path in bedGraph_files:
    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) != 6:
                continue  # Skip malformed lines

            chrom, start, end, methylation, meth_reads, unmeth_reads = parts

            # Normalize chromosome format
            if (
                chrom == snakemake.params["chromosome"]
                or chrom == "chr" + snakemake.params["chromosome"]
            ):
                try:
                    methylation = float(methylation)
                    meth_reads = int(meth_reads)
                    unmeth_reads = int(unmeth_reads)
                except ValueError:
                    continue  # Skip lines with conversion errors

                coverage = meth_reads + unmeth_reads
                if coverage == 0:
                    continue

                key = (chrom, start, end)
                if key not in bedGraph_entry:
                    bedGraph_entry[key] = [[methylation], meth_reads, unmeth_reads]
                else:
                    bedGraph_entry[key][0].append(methylation)
                    bedGraph_entry[key][1] += meth_reads
                    bedGraph_entry[key][2] += unmeth_reads

# Filter entries with at least 12 methylation values
bedGraph_entry = {
    key: value for key, value in bedGraph_entry.items() if len(value[0]) >= 12
}

# Filter entries with methylation range ≤ 20
bedGraph_entry = {
    key: value
    for key, value in bedGraph_entry.items()
    if max(value[0]) - min(value[0]) <= 20
}

# Write output to file
with open(snakemake.output[0], "w") as outfile:
    for key, values in bedGraph_entry.items():
        chrom, start, end = key
        meth_values, meth_reads, unmeth_reads = values
        total_reads = meth_reads + unmeth_reads

        avg_methylation = (meth_reads / total_reads) * 100 if total_reads != 0 else 0

        outfile.write(
            f"{chrom}\t{start}\t{end}\t{avg_methylation:.4f}\t"
            f"{meth_reads}\t{unmeth_reads}\n"
        )
