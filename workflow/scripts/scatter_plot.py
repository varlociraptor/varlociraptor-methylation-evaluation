import altair as alt
import pandas as pd
import numpy as np

import heapq


ref_dict = {}
vcf_dict = {}
true_dict = {}

ref_method = snakemake.wildcards.platform
# Create dictionaries for Bedgraph and VCF
# with open(snakemake.input["bedGraph"], 'r') as ref_file, open(snakemake.input["true_meth"][0], 'r') as truth_file:
with open(snakemake.input["calls"], 'r') as vcf_file, open(snakemake.input["bedGraph"], 'r') as ref_file, open(snakemake.input["true_meth"][0], 'r') as truth_file:

    match ref_method:
        case "Illumina_pe" | "Illumina_se":
            for line in ref_file:
                if not line.startswith("track"):
                    parts = line.strip().split('\t')
                    chrom, position, methylation_value = parts[0], (int(
                        parts[1]) + int(parts[2])) // 2, float(parts[3])
                    ref_dict[(chrom, position)] = methylation_value
        case "Nanopore":
            for line in ref_file:
                parts = line.strip().split()
                chrom, position, methylation_value = parts[0].removeprefix(
                    'chr'), int(parts[2]), float(parts[10])
                ref_dict[(chrom, position)] = methylation_value
        case "PacBio":
            for line in ref_file:
                if not line.startswith("track"):
                    parts = line.strip().split('\t')
                    chrom, position, methylation_value = parts[0], int(
                        parts[2]), float(parts[3])
                    ref_dict[(chrom, position)] = methylation_value
        case _:
            print("No valid method: ", ref_method)

    for line in vcf_file:
        if not line.startswith('#'):
            parts = line.strip().split('\t')
            chrom, pos = parts[0], int(parts[1])
            info_field = parts[8]
            af_index = info_field.split(":").index("AF")
            format_values = parts[9].split(":")
            af_value = float(format_values[af_index])
            vcf_dict[(chrom, pos)] = af_value

    for line in truth_file:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value = parts[0].replace(
                "chr", ""), (int(parts[1]) + int(parts[2])) // 2, float(parts[3])
            true_dict[(chrom, position)] = methylation_value


print("Ref_method data points: ", len(ref_dict))
print("Varlo data points: ", len(vcf_dict))
print("True data points: ", len(true_dict))
# ref_dict = {key: value for key, value in ref_dict.items() if value > 30}
# vcf_dict = {key: value for key, value in vcf_dict.items() if value > 0.3}


bedgraph_positions = [
    key for key in ref_dict if key in true_dict and key in vcf_dict]
bedgraph_meth_values = [ref_dict[key] for key in bedgraph_positions]

vcf_positions = [
    key for key in vcf_dict if key in ref_dict and key in true_dict]
vcf_af_values = [vcf_dict[key] * 100 for key in vcf_positions]


true_dict = dict(sorted(true_dict.items(), key=lambda item: item[0][1]))
true_positions = [
    key for key in true_dict if key in ref_dict and key in vcf_dict]
true_meth_values = [true_dict[key] for key in true_positions]


missing_positions1 = [
    key for key in ref_dict if key not in true_dict and key not in vcf_dict]
missing_positions2 = [
    key for key in vcf_dict if key not in ref_dict and key not in true_dict]
missing_positions3 = [
    key for key in true_dict if key not in ref_dict and key not in vcf_dict]


# with open(snakemake.output["test"], "w") as datei:
#     for i, el in enumerate(true_meth_values):
#         if abs(el - vcf_af_values[i]) > 20:
#             datei.write(str(true_positions[i]) + "\t" + str(el) + "\n")
#             datei.write("Bedgraph_value:" + str(bedgraph_meth_values[i]) + "\t Callsvcf:  " + str(vcf_af_values[i]) + "\n\n")
def euclidian_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
############################################################################################################################################################


line = alt.Chart(pd.DataFrame({'x': [1, 100], 'y': [1, 100]})).mark_line(color='red').encode(
    x='x:Q',
    y='y:Q'
)

# Die Datenpunkte
x_values = np.linspace(0, 100, num=101)
y_values = np.linspace(0, 100, num=101)


# Plot TrueMeth vs Varlociraptor
distances = [euclidian_distance(x, y, x, x)
             for x, y in zip(true_meth_values, vcf_af_values)]
sorted_list = sorted(distances, reverse=True)
# print("Max dist: ", sorted_list[:4])
deviation = sum(distances) / len(true_meth_values)


# Debug:
# for i, (x, y) in enumerate(zip(true_meth_values, vcf_af_values)):
#     if euclidian_distance(x, y, x, x) > 40 and x > 90:
#         print(euclidian_distance(x, y, x, x), x, y)
#         print(vcf_positions[i])

# Extra for plotting distances
rounded_distances = [int(round(d)) for d in distances]
data = pd.DataFrame({'Rounded_Distances': rounded_distances})
chart = alt.Chart(data).mark_bar().encode(
    x=alt.X('Rounded_Distances:O', axis=alt.Axis(title='distance')),
    y=alt.Y('count():Q', axis=alt.Axis(title='number'))
).properties(
    title='Distance distribution Varlociraptor'
)
chart.save(snakemake.output["dist_tv"], scale_factor=2.0)


data = pd.DataFrame({
    'TrueMeth': true_meth_values,
    'Varlociraptor': vcf_af_values
})

# Create a scatter plot with points colored by density
scatter = alt.Chart(data).mark_circle(opacity=0.2).encode(
    x='TrueMeth',
    y='Varlociraptor'
)


final_chart = (scatter + line).properties(
    width=400,
    height=400,
    title=f'TrueMeth vs. Varlociraptor (Deviation: {deviation})'
)
final_chart.save(snakemake.output["tv"], scale_factor=2.0)


# Plot RefMethod vs Varlociraptor
distances = [euclidian_distance(x, y, x, x)
             for x, y in zip(bedgraph_meth_values, vcf_af_values)]
sorted_list = sorted(distances, reverse=True)
deviation = sum(distances) / len(bedgraph_meth_values)

# Runde die Kommazahlen auf ganze Zahlen
rounded_distances = [int(round(d)) for d in distances]

# Erzeuge ein Pandas DataFrame
data = pd.DataFrame({'Rounded_Distances': rounded_distances})

# Erzeuge das Altair-Plot
chart = alt.Chart(data).mark_bar().encode(
    x=alt.X('Rounded_Distances:O', axis=alt.Axis(title='distance')),
    y=alt.Y('count():Q', axis=alt.Axis(title='number'))
).properties(
    title='Distance distribution RefMethod'
)

chart.save(snakemake.output["dist_rv"], scale_factor=2.0)


for i, (x, y) in enumerate(zip(bedgraph_meth_values, vcf_af_values)):
    if euclidian_distance(x, y, x, x) > 50:
        print(euclidian_distance(x, y, x, x), x, y)
        print(vcf_positions[i])


data = pd.DataFrame({
    'RefMethod': bedgraph_meth_values,
    'Varlociraptor': vcf_af_values
})

scatter = alt.Chart(data).mark_circle(opacity=0.2).encode(
    x='RefMethod',
    y='Varlociraptor'
)

final_chart = (scatter + line).properties(
    width=400,
    height=400,
    title=f'RefMethod vs. Varlociraptor (Deviation: {deviation})'
)
final_chart.save(snakemake.output["rv"], scale_factor=2.0)


# Plot TrueMeth vs RefMethod
distances = [euclidian_distance(x, y, x, x) for x, y in zip(
    bedgraph_meth_values, true_meth_values)]
sorted_list = sorted(distances, reverse=True)
deviation = sum(distances) / len(bedgraph_meth_values)
# Runde die Kommazahlen auf ganze Zahlen
rounded_distances = [int(round(d)) for d in distances]

# Erzeuge ein Pandas DataFrame
data = pd.DataFrame({'Rounded_Distances': rounded_distances})

# Erzeuge das Altair-Plot
chart = alt.Chart(data).mark_bar().encode(
    x=alt.X('Rounded_Distances:O', axis=alt.Axis(title='distance')),
    y=alt.Y('count():Q', axis=alt.Axis(title='number'))
).properties(
    title='Distance distribution RefMethod'
)

chart.save(snakemake.output["dist_rt"], scale_factor=2.0)


data = pd.DataFrame({
    'TrueMeth': true_meth_values,
    'RefMethod': bedgraph_meth_values
})

scatter = alt.Chart(data).mark_circle(opacity=0.2).encode(
    x='TrueMeth',
    y='RefMethod'
)

final_chart = (scatter + line).properties(
    width=400,
    height=400,
    title=f'TrueMeth vs. RefMethod (Deviation: {deviation})'
)
final_chart.save(snakemake.output["rt"], scale_factor=2.0)
