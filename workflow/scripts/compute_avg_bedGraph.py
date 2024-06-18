import pickle
import numpy as np


bedGraph_files = [snakemake.input[i] for i in range(len(snakemake.input))]

bedGraph_entry = {}


all_keys = set()


for file_path in bedGraph_files:
    file_keys = set()

    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            chrom, start, end, methylation, meth_reads, unmeth_reads = parts
            print("Chromosomes: ", chrom, snakemake.params["chromosome"])
            if chrom == snakemake.params["chromosome"] or chrom == "chr" + snakemake.params["chromosome"]:
                methylation = float(methylation)
                meth_reads = int(meth_reads)
                unmeth_reads = int(unmeth_reads)
                coverage = meth_reads + unmeth_reads
                if coverage > 0:
                    key = (chrom, start, end)

                    file_keys.add(key)

                    if key not in bedGraph_entry:
                        bedGraph_entry[key] = [
                            [methylation], meth_reads, unmeth_reads]
                    else:
                        bedGraph_entry[key][0].append(methylation)
                        bedGraph_entry[key][1] += meth_reads
                        bedGraph_entry[key][2] += unmeth_reads

#     if not all_keys:
#         all_keys = file_keys
#     else:
#         all_keys.intersection_update(file_keys)

# keys_to_remove = set(bedGraph_entry.keys()) - all_keys
# for key in keys_to_remove:
#     del bedGraph_entry[key]

bedGraph_entry = {key: value for key,
                  value in bedGraph_entry.items() if len(value[0]) >= 12}


bedGraph_entry = {
    key: value for key, value in bedGraph_entry.items()
    if max(value[0]) - min(value[0]) <= 20
}


# with open("bedGraph_entry.pkl", "wb") as outfile:
#     pickle.dump(bedGraph_entry, outfile)


with open(snakemake.output[0], "w") as outfile:
    for key, values in bedGraph_entry.items():
        chrom, start, end = key
        meth, meth_reads, unmeth_reads = values[0], values[1], values[2]
        # if np.std(meth) / (np.mean(meth) + 0.01) > 0.2:
        avg_methylation = (meth_reads / (meth_reads + unmeth_reads)) * \
            100 if meth_reads + unmeth_reads != 0 else 0
        outfile.write(f"{chrom}\t{start}\t{end}\t{avg_methylation:.4f}\t{
                      meth_reads:.0f}\t{unmeth_reads:.0f}\n")
