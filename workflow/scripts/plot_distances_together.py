import altair as alt
import pandas as pd
import numpy as np
import re
import os
from collections import Counter
import matplotlib.pyplot as plt


def compute_rmse(predictions, targets):
    squared_errors = [(x - y) ** 2 for x, y in zip(predictions, targets)]
    mean_squared_error = sum(squared_errors) / len(predictions)
    return np.sqrt(mean_squared_error)


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


def get_euclidian_distance(x1, y1, x2, y2):
    """
    Calculates the Euclidean distance between two points.

    Args:
        x1: x-coordinate of the first point.
        y1: y-coordinate of the first point.
        x2: x-coordinate of the second point.
        y2: y-coordinate of the second point.

    Returns:
        The Euclidean distance between the two points.
    """
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def plot_distances(tool1_name, tool2_name, true_meth, varlo_values, ref_values, true_values, bias_vals, output):
    """
    Plots the distribution of distances between predicted methylation values and truth methylation values for two tools.

    Args:
        tool1_name: Name of the first tool.
        tool2_name: Name of the second tool.
        true_meth: Name of the truth methylation column.
        tool1_values: List of methylation values predicted by the first tool.
        tool2_values: List of methylation values predicted by the second tool.
        true_values: List of truth methylation values.
        output: Path to the output file.
    """
    # Calculate distances
    tool_values_no_bias = [v for v, bias in zip(
        varlo_values, bias_vals) if bias == "normal"]
    true_values_no_bias = [v for v, bias in zip(
        true_values, bias_vals) if bias == "normal"]
    ref_values_no_bias = [v for v, bias in zip(
        ref_values, bias_vals) if bias == "normal"]
    print("Len: ", len(tool_values_no_bias))
    varlo_distances = [int(round(get_euclidian_distance(x, y, x, x)))
                       for x, y in zip(tool_values_no_bias, true_values_no_bias)]
    ref_distances = [int(round(get_euclidian_distance(x, y, x, x)))
                     for x, y in zip(ref_values_no_bias, true_values_no_bias)]

    print(len(tool_values_no_bias))
    print(len(true_values_no_bias))
    print(len(ref_values_no_bias))

    rmse_varlo = round(compute_rmse(
        tool_values_no_bias, true_values_no_bias), 2)
    rmse_ref = round(compute_rmse(ref_values_no_bias, true_values_no_bias), 2)

    # varlo_short = [v for i, v in enumerate(
    #     tool_values_no_bias) if varlo_distances[i] < 30]
    # true_short = [v for i, v in enumerate(
    #     true_values_no_bias) if varlo_distances[i] < 30]

    # rmse_varlo = np.round(compute_rmse(
    #     varlo_short, true_short), 2)
    # rmse_ref = np.round(compute_rmse(
    #     ref_values_no_bias, true_values_no_bias), 2)

    # Combine data into a DataFrame
    data = pd.DataFrame({
        tool1_name: varlo_distances,
        tool2_name: ref_distances,
    })

    melted_data = pd.melt(data, value_vars=[tool1_name, tool2_name],
                          var_name='Tool', value_name='Distanz')
    # melted_data = melted_data[melted_data['Distanz'] <= 30]

    # Base chart with shared y-axis encoding
    base = alt.Chart(melted_data).encode(
        x=alt.X('Distanz:Q', title='Distance in %'),
        y=alt.Y('count():Q', title='Count'),
        color=alt.Color('Tool:N', legend=alt.Legend(title='Tool'))
    )

    # Line chart
    line = base.mark_line().properties(
        title=tool1_name + ' rmse: ' + str(rmse_varlo) +
        ' / ' + tool2_name + ' rmse: ' + str(rmse_ref)
    )

    # Point chart
    points = base.mark_point()

    # Combine line and point charts
    histogram = line + points

    # Speichern oder Anzeigen des Histogramms
    histogram.save(output, scale_factor=2.0)

    # Speichern oder Anzeigen des Histogramms

    # varlo_data = data['Varlociraptor'].value_counts(
    # ).sort_index().values
    # dackel_data = data['MethylDackel'].value_counts(
    # ).sort_index().values

    # varlo_data = varlo_data[:20]
    # dackel_data = dackel_data[:20]

    # plt.figure(figsize=(10, 6))
    # plt.plot(range(1, len(varlo_data) + 1), varlo_data, marker='o',
    #          linestyle='-', color='b', label='Varlociraptor')
    # plt.plot(range(1, len(dackel_data) + 1), dackel_data, marker='s',
    #          linestyle='--', color='r', label='MethylDackel')

    # plt.xlabel('Distance in %')
    # plt.ylabel('Counts')
    # plt.legend()
    # plt.savefig(output)


varlo_dict = {}
ref_dict = {}
coverage_bin_dicts = {}
true_dict = {}
ref_dict = {}


with open(snakemake.input["true_meth"], 'r') as truth_file, open(snakemake.input["tool"], 'r') as tool_file, open(snakemake.input["ref_tool"], 'r') as ref_file:

    file_name = os.path.splitext(os.path.basename(ref_file.name))[0]
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

    # if file_name != 'calls':
    #     coverage_bin_dicts[0] = {'ref_dict': ref_dict,
    #                              'bias_dict': {}, 'prob_present': {}}

    for line in truth_file:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value = parts[0].replace(
                "chr", ""), (int(parts[1]) + int(parts[2])) // 2, float(parts[3])
            true_dict[(chrom, position)] = methylation_value

# Extract methylation values and truth values
tool_dict = coverage_bin_dicts[0]['tool_dict']
bias_dict = coverage_bin_dicts[0]['bias_dict']
prob_dict = coverage_bin_dicts[0]['prob_present']
print(len(list(tool_dict.keys())), len(list(true_dict.keys())))


# Berechne die Schnittmenge der Schlüssel von ref_dict, true_dict und ref_dict
all_cpg_positions = list(set(tool_dict.keys()) & set(
    true_dict.keys()) & set(ref_dict.keys()))

print(len(all_cpg_positions))


bias_pos = [
    bias_dict[key] if key in tool_dict and key in true_dict else "not in true" if key in tool_dict else "not in tool" for key in all_cpg_positions]

prob_pos = [prob_dict[key]
            if key in prob_dict else 1 for key in all_cpg_positions]

print("Bias: ", bias_pos.count("normal"))

tool_meth_values = [tool_dict.get(key, 0) for key in all_cpg_positions]
ref_meth_values = [ref_dict.get(key, 0) for key in all_cpg_positions]
true_meth_values = [true_dict.get(key, 0) for key in all_cpg_positions]


# Plot TrueMeth vs toolMethod
png_output = snakemake.output['png']


# Plot the distribution of distances
plot_distances("Varlociraptor", "Modkit", "TruthMeth",
               tool_meth_values, ref_meth_values, true_meth_values, bias_pos, png_output)
