from collections import defaultdict

combined_bedgraph_meth = defaultdict(float)
combined_bedgraph_unmeth = defaultdict(float)


bedgraph_files = snakemake.input

for bedgraph_file in bedgraph_files:
    with open(bedgraph_file, 'r') as bedgraph:
        for line in bedgraph:
            if not line.startswith("track"):
                parts = line.strip().split('\t')
                chrom, start, end, meth, unmeth = parts[0], int(
                    parts[1]), int(parts[2]), float(parts[4]), float(parts[5])
                combined_bedgraph_meth[(chrom, start, end)] += meth
                combined_bedgraph_unmeth[(chrom, start, end)] += unmeth

with open(snakemake.output[0], "w") as outfile:
    for (chrom, start, end), meth in dict(sorted(combined_bedgraph_meth.items(), key=lambda item: item[0][1])).items():
        unmeth = combined_bedgraph_unmeth[(chrom, start, end)]
        avg_methylation = (meth / (meth + unmeth)) * \
            100 if meth + unmeth != 0 else 0
        outfile.write(f"{chrom}\t{start}\t{end}\t{
                      avg_methylation:.4f}\t{meth:.0f}\t{unmeth:.0f}\n")
