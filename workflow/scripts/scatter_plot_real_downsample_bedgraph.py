import altair as alt
import pandas as pd

import heapq


dackel_dict = {}
vcf_dict = {}
true_dict = {}

# Create dictionaries for Bedgraph and VCF
with open(snakemake.input["calls"], 'r') as vcf_file, open(snakemake.input["bedGraph"], 'r') as dackel_file, open(snakemake.input["true_meth"], 'r') as truth_file:

    for line in dackel_file:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value = parts[0], (int(
                parts[1]) + int(parts[2])) // 2, float(parts[3])
            dackel_dict[(chrom, position)] = methylation_value

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


# dackel_dict = {key: value for key, value in dackel_dict.items() if value > 30}
# vcf_dict = {key: value for key, value in vcf_dict.items() if value > 0.3}


# bedgraph_positions = [key for key in dackel_dict if key in vcf_dict and key in true_dict]
bedgraph_positions = [key for key in dackel_dict if key in true_dict]
bedgraph_meth_values = [dackel_dict[key] for key in bedgraph_positions]
print(len(bedgraph_positions))

vcf_positions = [
    key for key in vcf_dict if key in dackel_dict and key in true_dict]
vcf_af_values = [vcf_dict[key] * 100 for key in vcf_positions]


true_dict = dict(sorted(true_dict.items(), key=lambda item: item[0][1]))
# true_positions = [key for key in true_dict if key in dackel_dict and key in vcf_dict]
true_positions = [key for key in true_dict if key in dackel_dict]
true_meth_values = [true_dict[key] for key in true_positions]
print(len(true_positions))


missing_positions1 = [
    key for key in dackel_dict if key not in vcf_dict and key not in true_dict]
missing_positions2 = [
    key for key in vcf_dict if key not in dackel_dict and key not in true_dict]
missing_positions3 = [
    key for key in true_dict if key not in dackel_dict and key not in vcf_dict]


# with open(snakemake.output["test"], "w") as datei:
#     for i, el in enumerate(true_meth_values):
#         if abs(el - vcf_af_values[i]) > 20:
#             datei.write(str(true_positions[i]) + "\t" + str(el) + "\n")
#             datei.write("Bedgraph_value:" + str(bedgraph_meth_values[i]) + "\t Callsvcf:  " + str(vcf_af_values[i]) + "\n\n")

############################################################################################################################################################


line = alt.Chart(pd.DataFrame({'x': [1, 100], 'y': [1, 100]})).mark_line(color='red').encode(
    x='x:Q',
    y='y:Q'
)

# Plot TrueMeth vs Varlociraptor
# deviation = sum(abs(x - y) for x, y in zip(true_meth_values, vcf_af_values))

# data = pd.DataFrame({
#     'TrueMeth': true_meth_values,
#     'Varlociraptor': vcf_af_values
# })

# scatter = alt.Chart(data).mark_circle(opacity=0.5).encode(
#     x='TrueMeth:Q',
#     y='Varlociraptor:Q'
# )

# final_chart = (scatter + line).properties(
#     width=400,
#     height=400,
#     title=f'TrueMeth vs. Varlociraptor (Deviation: {deviation})'
# )
# final_chart.save(snakemake.output["tv"], scale_factor=2.0)


# Plot MethylDackel vs Varlociraptor
# deviation = sum(abs(x - y) for x, y in zip(bedgraph_meth_values, vcf_af_values))


# data = pd.DataFrame({
#     'MethylDackel': bedgraph_meth_values,
#     'Varlociraptor': vcf_af_values
# })

# scatter = alt.Chart(data).mark_circle(opacity=0.5).encode(
#     x='MethylDackel:Q',
#     y='Varlociraptor:Q'
# )

# final_chart = (scatter + line).properties(
#     width=400,
#     height=400,
#     title=f'MethylDackel vs. Varlociraptor (Deviation: {deviation})'
# )
# final_chart.save(snakemake.output["dv"], scale_factor=2.0)


# Plot TrueMeth vs MethylDackel
deviation = (sum(abs(x - y) for x, y in zip(true_meth_values,
             bedgraph_meth_values))) / len(true_meth_values)

data = pd.DataFrame({
    'TrueMeth': true_meth_values,
    'MethylDackel': bedgraph_meth_values
})

scatter = alt.Chart(data).mark_circle(opacity=0.5).encode(
    x='TrueMeth:Q',
    y='MethylDackel:Q'
)

final_chart = (scatter + line).properties(
    width=400,
    height=400,
    title=f'TrueMeth vs. MethylDackel (Deviation: {deviation})'
)
final_chart.save(snakemake.output["td"], scale_factor=2.0)
