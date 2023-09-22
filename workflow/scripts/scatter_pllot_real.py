import altair as alt
import pandas as pd

import heapq



bedgraph_dict = {}
vcf_dict = {}
ch3_dict = {}

# Create dictionaries for Bedgraph and VCF
with open(snakemake.input["calls"], 'r') as vcf_file, open(snakemake.input["true_meth"], 'r') as bedgraph_file:

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
bedgraph_meth_values = [bedgraph_dict[key] for key in bedgraph_positions]

vcf_positions = [key for key in vcf_dict if key in bedgraph_dict]
vcf_af_values = [vcf_dict[key] * 100 for key in vcf_positions]


missing_positions = [key for key in bedgraph_dict if key not in vcf_dict]

for i in range (0, 7):
    println(len(bedgraph_meth_values), len(vcf_af_values))
    println(missing_positions)


# deviations = [abs(x - y) for x, y in zip(ch3_meth_values, vcf_af_values)]
# top_10_positions_with_deviation = heapq.nlargest(10, zip(deviations, ch3_positions))
# top_10_deviations, top_10_positions = zip(*top_10_positions_with_deviation)

# for i, (deviation, position) in enumerate(zip(top_10_deviations, top_10_positions), 1):
#     print(f"Position {i}: {position}, Deviation: {deviation}")




line = alt.Chart(pd.DataFrame({'x': [1, 100], 'y': [1, 100]})).mark_line(color='red').encode(
    x='x:Q',
    y='y:Q'
)




# Plot TrueMeth vs Varlociraptor
deviation = sum(abs(x - y) for x, y in zip(bedgraph_meth_values, vcf_af_values))


data = pd.DataFrame({
    'TrueMeth': bedgraph_meth_values,
    'Varlociraptor': vcf_af_values
})

scatter = alt.Chart(data).mark_circle(opacity=0.5).encode(
    x='TrueMeth:Q',
    y='Varlociraptor:Q'
)

final_chart = (scatter + line).properties(
    width=400,
    height=400, 
    title=f'TrueMeth vs. Varlociraptor (Deviation: {deviation})'
)
final_chart.save(snakemake.output["dv"], scale_factor=2.0) 

