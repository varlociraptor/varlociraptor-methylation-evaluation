import matplotlib.pyplot as plt


bedgraph_dict = {}
vcf_dict = {}

# Create dictionaries for Bedgraph and VCF
with open(snakemake.input["calls"], 'r') as vcf_file, open(snakemake.input["bedGraph"], 'r') as bedgraph_file:

    for line in bedgraph_file:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value = parts[0], (int(parts[1]) + int(parts[2])) // 2, int(parts[3])
            bedgraph_dict[(chrom, position)] = methylation_value

    for line in vcf_file:
        if not line.startswith('#'):
            parts = line.strip().split('\t')
            chrom, pos = parts[0], int(parts[1])
            info_field = parts[8]
            af_index = info_field.split(":").index("AF")
            format_values = parts[9].split(":")
            af_value = float(format_values[af_index])
            vcf_dict[(chrom, pos)] = af_value
 


bedgraph_positions = [key for key in bedgraph_dict if key in vcf_dict]
bedgraph_methylation_values = [bedgraph_dict[key] for key in bedgraph_positions]

vcf_positions = [key for key in vcf_dict if key in bedgraph_dict]
vcf_af_values = [vcf_dict[key] * 100 for key in vcf_positions]

missing_positions1 = [key for key in bedgraph_dict if key not in vcf_dict]
missing_positions2 = [key for key in vcf_dict if key not in bedgraph_dict]



# print(bedgraph_methylation_values)
# Create a scatterplot
fig, ax = plt.subplots()
ax.annotate(
    f"Missing positions in Varlociraptor:{missing_positions1}. \nMissing positions in MethylDackel:{missing_positions2}.",
    xy=(0.5, -0.1),  # Koordinaten der Textposition (x, y)
    xycoords='axes fraction',  # Verwenden Sie die Koordinatenachse des Plots
    fontsize=6,
    ha="center",
    va="center"
)
plt.plot([1, 100], [1, 100], label='Proportionale Ausgleichsgerade', color='red')

plt.scatter(vcf_af_values, bedgraph_methylation_values, alpha=1)
plt.xlabel('AF Values')
plt.ylabel('Methylation Values')
plt.title('MethylDackel vs. Varlociraptor')



plt.savefig(snakemake.output[0])
