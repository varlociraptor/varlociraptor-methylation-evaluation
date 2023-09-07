import matplotlib.pyplot as plt

# Pfad zur Bedgraph-Datei
bedgraph_path = snakemake.input["bedGraph"]
vcf_path = snakemake.input["calls"]

# Initialize dictionaries to store position-value mappings
bedgraph_dict = {}
vcf_dict = {}

# 2. Create dictionaries for Bedgraph and VCF
with open(vcf_path, 'r') as vcf_file, open(bedgraph_path, 'r') as bedgraph_file:

    for line in bedgraph_file:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value = parts[0], int(parts[1]) + int(parts[2]) // 2, int(parts[3])
            # Store the position and methylation value in the dictionary
            bedgraph_dict[(chrom, position)] = methylation_value

    for line in vcf_file:
        if not line.startswith('#'):
            parts = line.strip().split('\t')
            chrom, pos = parts[0], int(parts[1])
            info_field = parts[8]
            af_index = info_field.split(":").index("AF")
            format_values = parts[9].split(":")
            af_value = float(format_values[af_index])
            # Füge den AF-Wert zur Liste hinzu
            vcf_dict[(chrom, pos)] = af_value
 

""" # Print the dictionaries
print("Bedgraph Dictionary:")
for key, value in bedgraph_dict.items():
    print(f'Position: {key}, Methylation Value: {value}') """

# print("VCF Dictionary:")
# for key, value in vcf_dict.items():
#     print(f'Position: {key}, AF Value: {value}')



# 3. Scatterplot and remove redundant positions
bedgraph_positions = []
bedgraph_methylation_values = []
vcf_positions = []
vcf_af_values = []
test = True


for key, value in bedgraph_dict.items():
    if key in list(vcf_dict.keys()):
        bedgraph_positions.append(key)
        bedgraph_methylation_values.append(value)

for key, value in vcf_dict.items():
    if key in list(bedgraph_dict.keys()):

        vcf_positions.append(key)
        vcf_af_values.append(value * 100)


# print(bedgraph_methylation_values)
# Create a scatterplot
plt.scatter(vcf_af_values, bedgraph_methylation_values, alpha=0.5)
plt.xlabel('AF Values (Percentage)')
plt.ylabel('Methylation Values')
plt.title('Scatterplot: Methylation vs. AF')
plt.savefig(snakemake.output[0])
plt.show()