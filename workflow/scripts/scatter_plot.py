import altair as alt
import pandas as pd

import heapq



bedgraph_dict = {}
vcf_dict = {}
ch3_dict = {}

# Create dictionaries for Bedgraph and VCF
with open(snakemake.input["calls"], 'r') as vcf_file, open(snakemake.input["bedGraph"], 'r') as bedgraph_file, open(snakemake.input["true_meth"], 'r') as ch3_file:

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
 
    true_meth_data = ch3_file.readlines()
    for i in range(len(true_meth_data) - 1):
        parts1 = true_meth_data[i].strip().split('\t')
        parts2 = true_meth_data[i + 1].strip().split('\t')
        
        if parts1[3] == 'CG' and parts2[3] == 'CG' and parts1[2] == "+":
            chrom = parts1[0]
            position1 = int(parts1[1])  
            position2 = int(parts2[1]) 
            methylation1 = float(parts1[6]) 
            methylation2 = float(parts2[6])  

            avg_methylation = (methylation1 + methylation2) / 2
            ch3_dict[(chrom, position2)] = avg_methylation
            i += 1


bedgraph_positions = [key for key in bedgraph_dict if key in vcf_dict and key in ch3_dict]
bedgraph_meth_values = [bedgraph_dict[key] for key in bedgraph_positions]

vcf_positions = [key for key in vcf_dict if key in bedgraph_dict and key in ch3_dict]
vcf_af_values = [vcf_dict[key] * 100 for key in vcf_positions]

ch3_positions = [key for key in ch3_dict if key in bedgraph_dict and key in vcf_dict]
ch3_meth_values = [ch3_dict[key] * 100 for key in ch3_positions]

missing_positions1 = [key for key in bedgraph_dict if key not in vcf_dict and key not in ch3_dict]
missing_positions2 = [key for key in vcf_dict if key not in bedgraph_dict and key not in ch3_dict]
missing_positions3 = [key for key in ch3_dict if key not in bedgraph_dict and key not in vcf_dict]


# Berechnen Sie die Abweichungen zwischen ch3_meth_values und vcf_af_values
deviations = [abs(x - y) for x, y in zip(ch3_meth_values, vcf_af_values)]

# Verwenden Sie heapq.nlargest, um die Positionen mit den größten 10 Abweichungswerten zu erhalten
top_10_positions_with_deviation = heapq.nsmallest(30, zip(deviations, ch3_positions))

# Trennen Sie die Positionen und Abweichungswerte
top_10_deviations, top_10_positions = zip(*top_10_positions_with_deviation)

# Drucken Sie die Positionen mit den größten 10 Abweichungswerten
for i, (deviation, position) in enumerate(zip(top_10_deviations, top_10_positions), 1):
    print(f"Position {i}: {position}, Abweichung: {deviation}")



# print(bedgraph_methylation_values)
# Create a scatterplot
# Create the proportional straight line
line = alt.Chart(pd.DataFrame({'x': [1, 100], 'y': [1, 100]})).mark_line(color='red').encode(
    x='x:Q',
    y='y:Q'
)

# Plot TrueMeth vs Varlociraptor
deviation = sum(abs(x - y) for x, y in zip(ch3_meth_values, vcf_af_values))

data = pd.DataFrame({
    'TrueMeth': ch3_meth_values,
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
final_chart.save(snakemake.output["tv"], scale_factor=2.0) 



# Plot MethylDackel vs Varlociraptor

deviation = sum(abs(x - y) for x, y in zip(bedgraph_meth_values, vcf_af_values))


data = pd.DataFrame({
    'MethylDackel': bedgraph_meth_values,
    'Varlociraptor': vcf_af_values
})

scatter = alt.Chart(data).mark_circle(opacity=0.5).encode(
    x='MethylDackel:Q',
    y='Varlociraptor:Q'
)

final_chart = (scatter + line).properties(
    width=400,
    height=400, 
    title=f'MethylDackel vs. Varlociraptor (Deviation: {deviation})'
)
final_chart.save(snakemake.output["dv"], scale_factor=2.0) 


# Plot TrueMeth vs MethylDackel
deviation = sum(abs(x - y) for x, y in zip(ch3_meth_values, bedgraph_meth_values))

print(len(bedgraph_meth_values))
print(len(bedgraph_dict))

data = pd.DataFrame({
    'TrueMeth': ch3_meth_values,
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


