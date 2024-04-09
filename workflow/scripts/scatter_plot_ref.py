import altair as alt
import pandas as pd
import numpy as np
import re
import os
import pickle
import heapq


def entry_is_bias(prob, format_values):
    probs = re.split("[=;]", prob)
    probs_values = [probs[i] for i in [1, 3, 5]]
    try:
        prob_values = [float(value) for value in probs_values]
    except:
        return True
    min_value = min(prob_values)
    min_index = prob_values.index(min_value)
    bias = format_values[6:12]
    no_bias = all(b == '.' for b in bias)
    if min_index == 0 and not no_bias:
        return True
    return False


def euclidian_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def plot_distances(tool_name, true_meth, tool_values, true_values, bias_vals, output):
    distances = [euclidian_distance(x, y, x, x)
                 for x, y, bias in zip(tool_values, true_values, bias_vals) if not bias]
    print(distances)
    #  for x, y in zip(tool_values, true_values)]
    # sorted_list = sorted(distances, reverse=True)
    deviation = sum(distances) / len(true_meth_values)
    rounded_distances = [int(round(d)) for d in distances]
    data = pd.DataFrame({'Rounded_Distances': rounded_distances})
    chart = alt.Chart(data).mark_bar().encode(
        x=alt.X('Rounded_Distances:O', axis=alt.Axis(title='distance')),
        y=alt.Y('count():Q', axis=alt.Axis(title='number'))
    ).properties(
        title=f'Distance {tool_name} vs. {true_meth} (Deviation: {deviation})'
    )
    chart.save(output, scale_factor=2.0)
    return deviation


def plot_meth_vals(tool_name, true_meth, tool_values, true_values, bias_vals, output, deviation):
    data = pd.DataFrame({
        tool_name: tool_values,
        true_meth: true_values,
        'Bias': bias_vals
    })
    print(data)

    scatter = alt.Chart(data).mark_circle(opacity=0.2).encode(
        x=tool_name,
        y=true_meth,
        color='Bias'
    )

    final_chart = (scatter + line).properties(
        width=400,
        height=400,
        title=f'{tool_name} vs. {true_meth} (Deviation: {deviation})'
    )
    final_chart.save(output, scale_factor=2.0)


# truemeth[0], Illumina_pe
tool_dict = {}
true_dict = {}
bias_dict = {}
filename = ""


# with open(snakemake.input["calls"], 'r') as vcf_file, open(snakemake.input["tool"], 'r') as tool_file, open(snakemake.input["true_meth"][0], 'r') as truth_file:
with open(snakemake.input["true_meth"][0], 'r') as truth_file, open(snakemake.input["tool"], 'r') as tool_file:
    file_name = os.path.splitext(os.path.basename(tool_file.name))[0]

    if file_name == 'calls':  # Varlociraptor
        for line in tool_file:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom, pos = parts[0], int(parts[1])
                probs_field = parts[7]
                info_field = parts[8]
                af_index = info_field.split(":").index("AF")
                format_values = parts[9].split(":")
                af_value = float(format_values[af_index])
                tool_dict[(chrom, pos)] = af_value * 100
                if entry_is_bias(probs_field, format_values):
                    bias_dict[(chrom, pos)] = True

                else:
                    bias_dict[(chrom, pos)] = False
    # file_name = os.path.splitext(os.path.basename(tool_file))[0]

    if file_name == 'methylDackel':
        for line in tool_file:
            if not line.startswith("track"):
                parts = line.strip().split('\t')
                chrom, position, methylation_value = parts[0], (int(
                    parts[1]) + int(parts[2])) // 2, float(parts[3])
                tool_dict[(chrom, position)] = methylation_value

    if file_name == 'bsMap':
        for line in tool_file:
            if not line.startswith("chr"):
                parts = line.strip().split('\t')
                chrom, position, methylation_value = parts[0], int(
                    parts[1]), float(parts[4]) * 100
                tool_dict[(chrom, position)] = methylation_value

    if file_name == 'bismark':
        for line in tool_file:
            if not line.startswith("track"):
                parts = line.strip().split('\t')
                chrom, position, methylation_value = parts[0], (int(
                    parts[1]) + int(parts[2])) // 2, float(parts[3])
                tool_dict[(chrom, position)] = methylation_value

    if file_name == 'bisSNP':
        for line in tool_file:
            if not line.startswith("track"):
                parts = line.strip().split('\t')
                chrom, position, methylation_value = parts[0], int(
                    parts[2]), float(parts[3])
                tool_dict[(chrom, position)] = methylation_value

    if file_name == 'modkit':
        for line in tool_file:
            parts = line.strip().split()
            chrom, position, methylation_value = parts[0].removeptoolix(
                'chr'), int(parts[2]), float(parts[10])
            tool_dict[(chrom, position)] = methylation_value

    if file_name == 'pb_CpG_tools':
        for line in tool_file:
            if not line.startswith("track"):
                parts = line.strip().split('\t')
                chrom, position, methylation_value = parts[0], int(
                    parts[2]), float(parts[3])
                tool_dict[(chrom, position)] = methylation_value

    for line in truth_file:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value = parts[0].replace(
                "chr", ""), (int(parts[1]) + int(parts[2])) // 2, float(parts[3])
            true_dict[(chrom, position)] = methylation_value


print("tool_method data points: ", len(tool_dict))
print("True data points: ", len(true_dict))


# # Compute positions that are included in all tool methods
# common_tool_keys = list(tool_dicts[0][2].keys())
# # common_tool_keys = sorted(common_tool_keys, key=lambda x: x[1])
# for _, _, tool_dict in tool_dicts[1:]:
#     common_tool_keys = [c for c in common_tool_keys if c in tool_dict]

# tool_set = set(tool_dicts[0][2].keys())
# common_set = set(common_tool_keys)

cpg_positions = [key for key in tool_dict if key in true_dict]
# We cant keep all values as this leads to memory error
cpg_shortened = cpg_positions[0:100000]

if file_name != 'calls':
    # with open(snakemake.input['bias_pos'], "rb") as f:
    #     bias_pos = list((f.read()))
    #     bias_set = set(bias_pos)
    #     cpg_set = set(cpg_shortened)
    #     bias_pos = [True if pos in bias_set else False for pos in cpg_set]
    bias_pos = [False for _ in cpg_shortened]
else:
    bias_pos = [bias_dict[key] for key in cpg_shortened]


tool_dict = dict(sorted(tool_dict.items(), key=lambda item: item[0][1]))
bias_dict = dict(sorted(bias_dict.items(), key=lambda item: item[0][1]))
true_dict = dict(sorted(true_dict.items(), key=lambda item: item[0][1]))

tool_meth_values = [tool_dict[key] for key in cpg_shortened]
true_meth_values = [true_dict[key] for key in cpg_shortened]

line = alt.Chart(pd.DataFrame({'x': [1, 100], 'y': [1, 100]})).mark_line(color='red').encode(
    x='x:Q',
    y='y:Q'
)


# Plot TrueMeth vs toolMethod
deviation = plot_distances(file_name, 'TrueMeth', tool_meth_values,
                           true_meth_values, bias_pos, snakemake.output['tool_dist'])
deviation = plot_meth_vals(file_name, 'TrueMeth', tool_meth_values,
                           true_meth_values, bias_pos, snakemake.output['tool'], deviation)

if file_name == 'calls':
    with open(snakemake.output['bias_pos'], 'w') as f:
        f.write(str(bias_dict))
