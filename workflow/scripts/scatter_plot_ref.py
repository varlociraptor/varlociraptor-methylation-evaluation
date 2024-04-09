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


def plot_distances(x_name, y_name, vals1, vals2, bias_vals, output):
    distances = [euclidian_distance(x, y, x, x)
                 #  for x, y, bias in zip(vals1, vals2, bias_vals) if not bias]
                 for x, y in zip(vals1, vals2)]
    # sorted_list = sorted(distances, reverse=True)
    deviation = sum(distances) / len(true_meth_values)
    rounded_distances = [int(round(d)) for d in distances]
    data = pd.DataFrame({'Rounded_Distances': rounded_distances})
    chart = alt.Chart(data).mark_bar().encode(
        x=alt.X('Rounded_Distances:O', axis=alt.Axis(title='distance')),
        y=alt.Y('count():Q', axis=alt.Axis(title='number'))
    ).properties(
        title=f'Distance {x_name} vs. {y_name} (Deviation: {deviation})'
    )
    chart.save(output, scale_factor=2.0)
    return deviation


def plot_meth_vals(x_name, y_name, vals1, vals2, bias_vals, output, deviation):
    data = pd.DataFrame({
        x_name: vals1,
        y_name: vals2,
    })
    # 'Bias': bias_vals

    scatter = alt.Chart(data).mark_circle(opacity=0.2).encode(
        x=x_name,
        y=y_name,
    )
    # color='Bias'

    final_chart = (scatter + line).properties(
        width=400,
        height=400,
        title=f'{x_name} vs. {y_name} (Deviation: {deviation})'
    )
    final_chart.save(output, scale_factor=2.0)


# truemeth[0], Illumina_pe
ref_dict = {}
true_dict = {}

# All values Varlociraptor marks as biased
bias_dict = {}

# with open(snakemake.input["calls"], 'r') as vcf_file, open(snakemake.input["ref"], 'r') as ref_file, open(snakemake.input["true_meth"][0], 'r') as truth_file:
with open(snakemake.input["true_meth"][0], 'r') as truth_file, open(snakemake.input["ref"], 'r') as ref_file:
    file_name = os.path.splitext(os.path.basename(ref_file.name))[0]

    # file_name = os.path.splitext(os.path.basename(ref_file))[0]
    if file_name == 'methylDackel':
        for line in ref_file:
            if not line.startswith("track"):
                parts = line.strip().split('\t')
                chrom, position, methylation_value = parts[0], (int(
                    parts[1]) + int(parts[2])) // 2, float(parts[3])
                ref_dict[(chrom, position)] = methylation_value
    if file_name == 'bsMap':
        for line in ref_file:
            if not line.startswith("chr"):
                parts = line.strip().split('\t')
                chrom, position, methylation_value = parts[0], int(
                    parts[1]), float(parts[4]) * 100
                ref_dict[(chrom, position)] = methylation_value
    if file_name == 'bismark':
        for line in ref_file:
            if not line.startswith("track"):
                parts = line.strip().split('\t')
                chrom, position, methylation_value = parts[0], (int(
                    parts[1]) + int(parts[2])) // 2, float(parts[3])
                ref_dict[(chrom, position)] = methylation_value
    if file_name == 'bisSNP':
        for line in ref_file:
            if not line.startswith("track"):
                parts = line.strip().split('\t')
                chrom, position, methylation_value = parts[0], int(
                    parts[2]), float(parts[3])
                ref_dict[(chrom, position)] = methylation_value
    if file_name == 'modkit':
        for line in ref_file:
            parts = line.strip().split()
            chrom, position, methylation_value = parts[0].removeprefix(
                'chr'), int(parts[2]), float(parts[10])
            ref_dict[(chrom, position)] = methylation_value
    if file_name == 'pb_CpG_tools':

        for line in ref_file:
            if not line.startswith("track"):
                parts = line.strip().split('\t')
                chrom, position, methylation_value = parts[0], int(
                    parts[2]), float(parts[3])
                ref_dict[(chrom, position)] = methylation_value

    for line in truth_file:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value = parts[0].replace(
                "chr", ""), (int(parts[1]) + int(parts[2])) // 2, float(parts[3])
            true_dict[(chrom, position)] = methylation_value


print("Ref_method data points: ", len(ref_dict))
print("True data points: ", len(true_dict))


# # Compute positions that are included in all ref methods
# common_ref_keys = list(ref_dicts[0][2].keys())
# # common_ref_keys = sorted(common_ref_keys, key=lambda x: x[1])
# for _, _, ref_dict in ref_dicts[1:]:
#     common_ref_keys = [c for c in common_ref_keys if c in ref_dict]

# ref_set = set(ref_dicts[0][2].keys())
# common_set = set(common_ref_keys)

cpg_positions = [
    key for key in ref_dict if key in true_dict]
ref_dict = dict(sorted(ref_dict.items(), key=lambda item: item[0][1]))

ref_meth_values = [ref_dict[key] for key in cpg_positions]

with open(snakemake.input['bias_pos'], "rb") as f:
    bias_pos = list((f.read()))
bias_set = set(bias_pos)
cpg_set = set(cpg_positions)
bias_dict = {(pos): (True if pos in bias_set
                     else False) for pos in cpg_set}

bias_dict = dict(sorted(bias_dict.items(), key=lambda item: item[0][1]))
true_dict = dict(sorted(true_dict.items(), key=lambda item: item[0][1]))
true_meth_values = [true_dict[key] for key in cpg_positions]


line = alt.Chart(pd.DataFrame({'x': [1, 100], 'y': [1, 100]})).mark_line(color='red').encode(
    x='x:Q',
    y='y:Q'
)


# Plot TrueMeth vs RefMethod
deviation = plot_distances(file_name, 'TrueMeth', ref_meth_values,
                           true_meth_values, bias_dict, snakemake.output['ref_dist'])
deviation = plot_meth_vals(file_name, 'TrueMeth', ref_meth_values,
                           true_meth_values, bias_dict, snakemake.output['ref'], deviation)
