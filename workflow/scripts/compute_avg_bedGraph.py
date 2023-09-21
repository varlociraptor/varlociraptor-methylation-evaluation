import pickle
import numpy as np


bedGraph_files = [snakemake.input[i] for i in range(len(snakemake.input))]

bedGrap_entry = {}


all_keys = set()


for file_path in bedGraph_files:
    file_keys = set()
    
    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            chrom, start, end, methylation, meth_reads, unmeth_reads = parts

            methylation = float(methylation)
            meth_reads = int(meth_reads)
            unmeth_reads = int(unmeth_reads)

            key = (chrom, start, end)
            
            file_keys.add(key)

            if key not in bedGrap_entry:
                bedGrap_entry[key] = [[methylation], meth_reads, unmeth_reads]
            else:
                bedGrap_entry[key][0].append(methylation)
                bedGrap_entry[key][1] += meth_reads
                bedGrap_entry[key][2] += unmeth_reads
            print(bedGrap_entry[key], file_path)


    if not all_keys:
        all_keys = file_keys
    else:
        all_keys.intersection_update(file_keys)

keys_to_remove = set(bedGrap_entry.keys()) - all_keys
for key in keys_to_remove:
    del bedGrap_entry[key]

# with open("bedGraph_entry.pkl", "wb") as outfile:
#     pickle.dump(bedGrap_entry, outfile)


with open(snakemake.output[0], "w") as outfile:
    for key, values in bedGrap_entry.items():
        chrom, start, end = key
        meth, meth_reads, unmeth_reads = values[0], values[1], values[2]
        if np.std(meth) / (np.mean(meth) + 0.01) > 0.2:
            avg_methylation = (meth_reads / (meth_reads + unmeth_reads)) * 100 if meth_reads + unmeth_reads != 0 else 0
            outfile.write(f"{chrom}\t{start}\t{end}\t{avg_methylation:.4f}\t{meth_reads:.0f}\t{unmeth_reads:.0f}\n")
