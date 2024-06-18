import altair as alt
import pandas as pd
import numpy as np
import re
import os
import pickle
import heapq
from vega_datasets import data
import matplotlib.pyplot as plt


def get_bin(coverage):
    return min(int(coverage / 10), snakemake.params['cov_bins'] - 1)


def get_prob_present(info):
    pattern = r"PROB_PRESENT=([0-9.]+)"
    match = re.search(pattern, info)
    try:
        prob_present_value = float(match.group(1))
    except:
        return 0  # No prob value given because of bias
    return 10 ** (-prob_present_value / 10)


def assign_coverage_bin(coverage):
    return (coverage - 1) // snakemake.params["coverage_bin"]


def compute_bias(prob, format_values):
    probs = re.split("[=;]", prob)
    probs_values = [probs[i] for i in [1, 3, 5]]
    try:
        prob_values = [float(value) for value in probs_values]
    except:
        return 'bias'
    min_value = min(prob_values)
    min_index = prob_values.index(min_value)

    bias_labels = ['SB', 'ROB', 'RPB', 'SCB', 'HE', 'ALB']

    if any(value != '.' for value in format_values[6:12]):
        bias = bias_labels[format_values[6:12].index('.')]
    else:
        bias = 'normal'

    if min_index != 0:
        bias = 'normal'
    return bias


def compute_rmse(predictions, targets):
    squared_errors = [(x - y) ** 2 for x, y in zip(predictions, targets)]
    mean_squared_error = sum(squared_errors) / len(predictions)
    return np.sqrt(mean_squared_error)


def euclidian_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def plot_distances(tool_name, true_meth, tool_values, true_values, bias_vals, output):
    tool_values_no_bias = [v for v, bias in zip(
        tool_values, bias_vals) if bias == "normal"]
    true_values_no_bias = [v for v, bias in zip(
        true_values, bias_vals) if bias == "normal"]
    distances = [euclidian_distance(x, y, x, x)
                 for x, y in zip(tool_values_no_bias, true_values_no_bias)]
    # #  for x, y in zip(tool_values, true_values)]
    # # sorted_list = sorted(distances, reverse=True)

    rmse = compute_rmse(tool_values_no_bias, true_values_no_bias)

    # #  Find positions with largest distances
    # positions_no_bias = [pos for pos, bias in zip(
    #     all_cpg_positions, bias_vals) if bias == "normal"]

    # top_n = 500  # Adjust this value to get more or fewer top positions
    # largest_distances = heapq.nlargest(
    #     top_n, zip(distances, positions_no_bias))

    # for dist, pos in largest_distances:
    #     if pos in true_dict.keys():
    #         ref_value = ref_dict.get(pos, "N/A")
    #         varlo_value = varlo_dict.get(pos, "N/A")
    #         true_value = true_dict.get(pos, "N/A")

    #         print("Position: ", pos, "Distance: ", dist,
    #               "\n\tRef: ", ref_value,
    #               "\tVarlo: ", varlo_value,
    #               "\tTrue: ", true_value)

    #         if abs(true_value - ref_value) < abs(true_value - varlo_value):
    #             print("\tRef")
    #         else:
    #             print("\tVarlo")
    # Compute distances

    # distances_varlo_true = [compute_rmse(varlo_dict.get(
    #     pos, 21), true_dict.get(pos, 21)) for pos in all_cpg_positions]
    # distances_ref_true = [compute_rmse(ref_dict.get(
    #     pos, 21), true_dict.get(pos, 21)) for pos in all_cpg_positions]

    # # Calculate difference between distances
    # distance_diffs = [(pos, d1 - d2) for pos, d1,
    #                   d2 in zip(all_cpg_positions, distances_varlo_true, distances_ref_true)]

    # # Sort by difference in distances
    # sorted_diffs = sorted(distance_diffs, key=lambda x: x[1], reverse=True)

    # # Print positions where varlo distance is significantly worse than ref distance
    # for pos, diff in sorted_diffs[0:10]:
    #     ref_value = ref_dict.get(pos, "N/A")
    #     varlo_value = varlo_dict.get(pos, "N/A")
    #     true_value = true_dict.get(pos, "N/A")
    #     dist_varlo_true = distances_varlo_true[all_cpg_positions.index(
    #         pos)]
    #     dist_ref_true = distances_ref_true[all_cpg_positions.index(pos)]

    #     print("Position: ", pos, "Distance Difference: ", diff,
    #           "\n\tRef: ", ref_value, "Distance: ", dist_ref_true,
    #           "\tVarlo: ", varlo_value, "Distance: ", dist_varlo_true,
    #           "\tTrue: ", true_value)

    print(0)

    distances_varlo_true = [(pos, abs(varlo_dict.get(
        pos, 0) - true_dict.get(pos, 0))) for pos in all_cpg_positions]
    print(1)
    distances_ref_true = [(pos,
                           abs(ref_dict.get(pos, 0) - true_dict.get(pos, 0))) for pos in all_cpg_positions]
    print(2)
    distances_diff = [(pos, x - y) for (pos, x), (pos, y)
                      in zip(distances_varlo_true, distances_ref_true)]
    sorted_diffs = sorted(distances_diff, key=lambda x: x[1], reverse=True)

    # Sortiere die Positionen basierend auf dem Unterschied zwischen den Distanzen
    # positions_sorted_by_diff = sorted(all_cpg_positions, key=lambda pos: distances_varlo_true[all_cpg_positions.index(
    # pos)] - distances_ref_true[all_cpg_positions.index(pos)], reverse=True)

    print(3)

    # Gib die Positionen aus, bei denen die Distanz zwischen varlo und true am größten ist im Gegensatz zur Distanz zwischen ref und true
    for (pos, dist) in sorted_diffs[:10]:  # Top 10 Positionen
        ref_value = ref_dict.get(pos, "N/A")
        varlo_value = varlo_dict.get(pos, "N/A")
        true_value = true_dict.get(pos, "N/A")

        # dist_varlo_true = abs(varlo_value - true_value)
        # dist_ref_true = abs(ref_value - true_value)
        print("Position: ", pos,  "Distance: ", dist,
              "\n\tRef: ", ref_value,
              "\tVarlo: ", varlo_value,
              "\tTrue: ", true_value)

        # if dist_varlo_true > dist_ref_true:
        #     print("\tVarlo schlechter als Ref")
        # else:
        #     print("\tRef schlechter als Varlo")

    rounded_distances = [int(round(d)) for d in distances]
    data = pd.DataFrame({'Rounded_Distances': rounded_distances})
    chart = alt.Chart(data).mark_bar().encode(
        x=alt.X('Rounded_Distances:O', axis=alt.Axis(title='distance')),
        y=alt.Y('count():Q', axis=alt.Axis(title='number'))
    ).properties(
        title=f'Distance {tool_name} vs. {true_meth} (rmse: {rmse})'
    )
    chart.save(output, scale_factor=2.0)
    return rmse


def plot_meth_vals(tool_name, true_meth, tool_values, true_values, bias_vals, prob_vals, all_cpg_positions, png_output, html_output, rmse):
    chromosome = [chrom for (chrom, _) in all_cpg_positions]
    position = [position for (_, position) in all_cpg_positions]
    data = pd.DataFrame({
        tool_name: tool_values,
        true_meth: true_values,
        'Bias': bias_vals,
        'Prob': prob_vals,
        "chromosome": chromosome,
        "position": position
    })

    filtered_data = data[data['Bias'] == 'normal']

    # # Hexbins
    # plt.hexbin(filtered_data[tool_name], filtered_data[true_meth],
    #            gridsize=20, cmap='viridis', bins='log', mincnt=-1)

    # plt.colorbar(label='log10(N)')
    # plt.xlabel(tool_name)
    # plt.ylabel(true_meth)
    # plt.savefig(png_output, dpi=300, bbox_inches='tight')
    # plt.savefig(html_output, format='svg', dpi=300, bbox_inches='tight')

    print(len(data))
    if len(data) >= 700000:
        data = data.sample(n=700000)

    line = alt.Chart(pd.DataFrame({'x': [1, 100], 'y': [1, 100]})).mark_line(color='red').encode(
        x='x:Q',
        y='y:Q'
    )

    scatter = alt.Chart(data).mark_circle().encode(
        x=tool_name,
        y=true_meth,
        color="Bias:N",
        opacity="Prob:Q",
        tooltip=['chromosome:N', 'position:N']
    ).interactive()

    final_chart = (scatter + line).properties(
        width=400,
        height=400,
        title=alt.Title(f'{tool_name} vs. {true_meth}',
                        subtitle=f'Rmse: {rmse}, Data_points: {bias_vals.count("normal")}')
    ).interactive()

    final_chart.save(png_output, scale_factor=2.0)
    final_chart.save(html_output, scale_factor=2.0)


print(-3)
alt.data_transformers.enable("vegafusion")
# truemeth[0], Illumina_pe
varlo_dict = {}
ref_dict = {}
bias_dict = {}
coverage_bin_dicts = {}
true_dict = {}
filename = ""


# with open(snakemake.input["calls"], 'r') as vcf_file, open(snakemake.input["tool"], 'r') as tool_file, open(snakemake.input["true_meth"][0], 'r') as truth_file:
# with open(snakemake.input["true_meth"][0], 'r') as truth_file, open(snakemake.input["tool"], 'r') as tool_file:
with open(snakemake.input["true_meth"], 'r') as truth_file, open(snakemake.input["tool"], 'r') as tool_file, open(snakemake.input["ref_tool"], 'r') as ref_file:
    file_name = os.path.splitext(os.path.basename(tool_file.name))[0]
    if file_name == 'calls':  # Varlociraptor
        for line in tool_file:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom, pos, info_field, format_field, values = parts[0], int(
                    parts[1]), parts[7], parts[8], parts[9].split(":")

                format_fields = format_field.split(":")
                dp_index = format_fields.index("DP")
                af_index = format_fields.index("AF")

                meth_rate = float(values[af_index])
                coverage = float(values[dp_index])

                cov_bin = get_bin(coverage)
                if cov_bin not in coverage_bin_dicts:
                    coverage_bin_dicts[cov_bin] = {
                        'tool_dict': {}, 'bias_dict': {}, 'prob_present': {}}

                chrom_pos = (chrom, pos)
                coverage_bin_dicts[cov_bin]['tool_dict'][chrom_pos] = meth_rate * 100
                coverage_bin_dicts[cov_bin]['prob_present'][chrom_pos] = get_prob_present(
                    info_field)
                coverage_bin_dicts[cov_bin]['bias_dict'][chrom_pos] = compute_bias(
                    info_field, values)

    for line in truth_file:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value = parts[0].replace(
                "chr", ""), (int(parts[1]) + int(parts[2])) // 2, float(parts[3])
            true_dict[(chrom, position)] = methylation_value

# Nur fuer debugging reasons, um pb zu plotten
    for line in ref_file:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value = parts[0], int(
                parts[2]), float(parts[3])
            ref_dict[(chrom, position)] = methylation_value

print(-1)
for cov_bin, bin_dicts in coverage_bin_dicts.items():

    varlo_dict = bin_dicts['tool_dict']
    bias_dict = bin_dicts['bias_dict']
    prob_dict = bin_dicts['prob_present']

    all_cpg_positions = list(set(varlo_dict.keys()) & set(
        ref_dict.keys()) & set(true_dict.keys()))
    # all_cpg_positions = list(all_cpg_positions)[:200000]
    # all_cpg_positions = list(all_cpg_positions)
    if file_name != 'calls':
        bias_pos = [
            "normal" if key in varlo_dict and key in ref_dict else "not in true" if key in varlo_dict else "not in tool" for key in all_cpg_positions]
        prob_pos = [1 for _ in all_cpg_positions]
    else:
        bias_pos = [
            bias_dict[key] if key in varlo_dict and key in ref_dict else "not in true" if key in varlo_dict else "not in tool" for key in all_cpg_positions]

        prob_pos = [prob_dict[key]
                    if key in prob_dict else 1 for key in all_cpg_positions]

    tool_meth_values = [varlo_dict.get(key, 0) for key in all_cpg_positions]
    true_meth_values = [ref_dict.get(key, 0) for key in all_cpg_positions]
    # Plot TrueMeth vs toolMethod
    png_output = snakemake.output['png']
    html_output = snakemake.output['html']
    dist_output = snakemake.output['tool_dist']

    # If we want to have multiple plots for different coverages we have to chose the correct one
    if isinstance(png_output, list):
        png_output = png_output[cov_bin]
        html_output = html_output[cov_bin]
        dist_output = dist_output[cov_bin]

    if file_name == 'calls':
        file_name = "Varlociraptor"
    # print("Test4", tool_meth_values, true_meth_values)

    rmse = plot_distances(file_name, 'TrueMeth', tool_meth_values,
                          true_meth_values, bias_pos, dist_output)
    rmse = plot_meth_vals(file_name, 'TrueMeth', tool_meth_values,
                          true_meth_values, bias_pos, prob_pos, all_cpg_positions, png_output, html_output, rmse)

    # if file_name == 'calls':
    #     with open(snakemake.output['bias_pos'], 'w') as f:
    #         f.write(str(bias_dict))
