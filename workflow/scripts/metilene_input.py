from collections import defaultdict

call_files = snakemake.input
platforms = [path.split('/')[1] for path in call_files]
meth_ratios = defaultdict(lambda: defaultdict(dict))


print("All files", call_files)
for i, file in enumerate(call_files):
    with open(file, 'r') as calls:

        for line in calls:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom, pos, info_field, format_field, values = parts[0], int(
                    parts[1]), parts[7], parts[8], parts[9].split(":")

                af_index = format_field.split(":").index("AF")
                methylation_rate = float(values[af_index])
                meth_ratios[(chrom, pos)][platforms[i]] = methylation_rate

meth_ratios = dict(
    sorted(meth_ratios.items(), key=lambda item: (item[0][0], item[0][1])))

# Schreibe Methylierungsraten in Datei
with open(str(snakemake.output), "w") as outfile:

    # Schreibe Header
    outfile.write("chr\tpos\t")
    outfile.write(
        "\t".join(platforms))
    outfile.write("\n")

    # Schreibe Methylierungsraten
    for (chrom, pos), _ in meth_ratios.items():
        outfile.write(f"{chrom}\t{pos}\t")
        for platform in platforms:
            methylation_rate = meth_ratios[(chrom, pos)].get(platform, '-')
            outfile.write(f"{methylation_rate}\t")
        outfile.write("\n")
