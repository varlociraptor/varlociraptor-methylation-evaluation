import altair as alt
import pandas as pd
import numpy as np
import re
import os
import csv
import pickle
import heapq
from vega_datasets import data
import matplotlib.pyplot as plt


def compute_precision_recall(file_name, cov_bin, data):
    # Define methylation threshold for a locus to be considered methylated
    methylation_threshold = 5
    # Convert methylation values to binary classifications (1 for methylated, 0 for unmethylated)
    tool_methylated = [
        1 if meth > methylation_threshold else 0 for meth in data[file_name]
    ]
    true_methylated = [
        1 if meth > methylation_threshold else 0 for meth in data["TrueMeth"]
    ]

    # Calculate True Positives (TP), False Positives (FP), and False Negatives (FN)
    TP = sum(1 for t, gt in zip(tool_methylated, true_methylated) if t == 1 and gt == 1)
    FP = sum(1 for t, gt in zip(tool_methylated, true_methylated) if t == 1 and gt == 0)
    FN = sum(1 for t, gt in zip(tool_methylated, true_methylated) if t == 0 and gt == 1)

    # Calculate Precision and Recall
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0

    # Calculate coverages (example based on your provided code, adjust as needed)
    # Assuming cov_bin and snakemake.params["cov_bin_size"] are defined in your environment
    cov_bin_size = snakemake.params["cov_bin_size"]

    coverages = (
        str(cov_bin * cov_bin_size)
        + "-"
        + str(cov_bin * cov_bin_size + cov_bin_size - 1)
    )

    return precision, recall, coverages


def save_precision_recall(tool, coverage, no_sites, precision, recall, filename):
    # Datei anlegen mit Header, wenn sie nicht existiert
    try:
        with open(filename, "x", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["tool", "coverage", "number_sites", "precision", "recall"])
    except FileExistsError:
        pass  # Datei existiert bereits, also Header nicht neu schreiben

    # Ergebnisse anhängen
    with open(filename, "a", newline="") as file:
        print(file)
        writer = csv.writer(file)
        writer.writerow([tool, coverage, no_sites, precision, recall])


def get_bin(coverage):
    return min(
        int(coverage / snakemake.params["cov_bin_size"]),
        snakemake.params["cov_bins"] - 1,
    )


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
        return "bias"
    min_value = min(prob_values)
    min_index = prob_values.index(min_value)

    bias_labels = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]

    if any(value != "." for value in format_values[6:12]):
        bias = bias_labels[format_values[6:12].index(".")]
    else:
        bias = "normal"

    if min_index != 0:
        bias = "normal"
    return bias


def fill_coverage_bin_dict(coverage, chrom_pos, methylation):
    cov_bin = get_bin(coverage)

    if cov_bin not in coverage_bin_dicts:
        coverage_bin_dicts[cov_bin] = {
            "tool_dict": {},
            "bias_dict": {},
        }
    coverage_bin_dicts[cov_bin]["tool_dict"][chrom_pos] = methylation


def compute_rmse(predictions, targets):
    squared_errors = [(x - y) ** 2 for x, y in zip(predictions, targets)]
    mean_squared_error = sum(squared_errors) / len(predictions)
    return np.sqrt(mean_squared_error)


def euclidian_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def plot_distances(tool_name, true_meth, tool_values, true_values, bias_vals, output):
    tool_values_no_bias = [
        v for v, bias in zip(tool_values, bias_vals) if bias == "normal"
    ]
    true_values_no_bias = [
        v for v, bias in zip(true_values, bias_vals) if bias == "normal"
    ]
    distances = [
        euclidian_distance(x, y, x, x)
        for x, y in zip(tool_values_no_bias, true_values_no_bias)
    ]

    rmse = compute_rmse(tool_values_no_bias, true_values_no_bias)

    rounded_distances = [int(round(d)) for d in distances]
    data = pd.DataFrame({"Rounded_Distances": rounded_distances})

    chart = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            x=alt.X("Rounded_Distances:O", axis=alt.Axis(title="distance")),
            y=alt.Y("count():Q", axis=alt.Axis(title="number")),
        )
        .properties(title=f"Distance {tool_name} vs. {true_meth} (rmse: {rmse})")
    )
    print(output)
    chart.save(output, scale_factor=2.0)
    return rmse


def plot_meth_vals(
    tool_name,
    true_meth,
    data,
    output,
    rmse,
):

    if len(data) >= 700000:
        data = data.sample(n=700000)

    line = (
        alt.Chart(pd.DataFrame({"x": [1, 100], "y": [1, 100]}))
        .mark_line(color="red")
        .encode(x="x:Q", y="y:Q")
    )

    scatter = (
        alt.Chart(data)
        .mark_circle()
        .encode(
            x=tool_name,
            y=true_meth,
            # color="Bias:N",
            # opacity="Prob:Q",
            opacity=alt.value(0.3),
            tooltip=["chromosome:N", "position:N"],
        )
        .interactive()
    )

    final_chart = (
        (scatter + line)
        .properties(
            width=400,
            height=400,
            title=alt.Title(
                f"{tool_name} vs. {true_meth}",
                subtitle=f"Rmse: {rmse}, datapoints: {len(data)}",
            ),
        )
        .interactive()
    )
    print(output)
    final_chart.save(output, scale_factor=2.0)


alt.data_transformers.enable("vegafusion")
# truemeth[0], Illumina_pe
tool_dict = {}
true_dict = {}
bias_dict = {}
coverage_bin_dicts = {}
filename = ""

# I need these to compute precision and recall (Precision: #True positives/ #All positives tool) (Recall: #True positives/ #All positives truth)
number_positives_tool = 0
number_positives_truth = 0

# with open(snakemake.input["calls"], 'r') as vcf_file, open(snakemake.input["tool"], 'r') as tool_file, open(snakemake.input["true_meth"][0], 'r') as truth_file:
with open(snakemake.input["true_meth"][0], "r") as truth_file, open(
    snakemake.input["tool"], "r"
) as tool_file:
    file_name = os.path.splitext(os.path.basename(tool_file.name))[0]
    if file_name == "calls":  # Varlociraptor
        for line in tool_file:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                chrom, pos, alternative, info_field, format_field, values = (
                    parts[0],
                    int(parts[1]),
                    parts[4],
                    parts[7],
                    parts[8],
                    parts[9].split(":"),
                )
                if alternative == "<METH>":
                    format_fields = format_field.split(":")
                    dp_index = format_fields.index("DP")
                    af_index = format_fields.index("AF")

                    meth_rate = float(values[af_index])
                    coverage = float(values[dp_index])

                    fill_coverage_bin_dict(coverage, (chrom, pos), meth_rate * 100)

                    cov_bin = get_bin(coverage)
                    coverage_bin_dicts[cov_bin]["bias_dict"][(chrom, pos)] = (
                        compute_bias(info_field, values)
                    )

    if file_name == "methylDackel":
        for line in tool_file:
            if not line.startswith("track"):
                parts = line.strip().split("\t")
                chrom, position, methylation_value, coverage = (
                    parts[0],
                    (int(parts[1]) + int(parts[2])) // 2,
                    float(parts[3]),
                    (int(parts[4]) + int(parts[5])),
                )
                fill_coverage_bin_dict(coverage, (chrom, position), methylation_value)

    if file_name == "bsMap":
        for line in tool_file:
            if not line.startswith("chr"):
                parts = line.strip().split("\t")
                chrom, position, methylation_value, coverage = (
                    parts[0],
                    int(parts[1]),
                    float(parts[4]) * 100,
                    float(parts[5]),
                )
                fill_coverage_bin_dict(coverage, (chrom, position), methylation_value)

    if file_name == "bismark":
        for line in tool_file:
            if not line.startswith("track"):
                parts = line.strip().split("\t")
                chrom, position, methylation_value = (
                    parts[0],
                    (int(parts[1]) + int(parts[2])) // 2,
                    float(parts[3]),
                )
                fill_coverage_bin_dict(coverage, (chrom, position), methylation_value)

    if file_name == "bisSNP":
        for line in tool_file:
            if not line.startswith("track"):
                parts = line.strip().split("\t")
                chrom, position, methylation_value = (
                    parts[0],
                    int(parts[2]),
                    float(parts[3]),
                )
                fill_coverage_bin_dict(coverage, (chrom, position), methylation_value)

    if file_name == "modkit":
        for line in tool_file:
            parts = line.strip().split()
            chrom, position, methylation_value, coverage = (
                parts[0].removeprefix("chr"),
                int(parts[2]),
                float(parts[10]),
                float(parts[10]),
            )
            fill_coverage_bin_dict(coverage, (chrom, position), methylation_value)

    if file_name == "pb_CpG_tools":
        for line in tool_file:
            if not line.startswith("track"):
                parts = line.strip().split("\t")
                chrom, position, methylation_value, coverage = (
                    parts[0],
                    int(parts[2]),
                    float(parts[3]),
                    int(parts[5]),
                )
                fill_coverage_bin_dict(coverage, (chrom, position), methylation_value)

    for line in truth_file:
        if not line.startswith("track"):
            parts = line.strip().split("\t")
            chrom, position, methylation_value = (
                parts[0].replace("chr", ""),
                (int(parts[1]) + int(parts[2])) // 2,
                float(parts[3]),
            )
            true_dict[(chrom, position)] = methylation_value


for cov_bin, bin_dicts in coverage_bin_dicts.items():
    tool_dict = bin_dicts["tool_dict"]
    bias_dict = bin_dicts["bias_dict"]

    all_cpg_positions = list(set(tool_dict.keys()) | set(true_dict.keys()))
    # all_cpg_positions = list(all_cpg_positions)[:200000]
    # all_cpg_positions = list(all_cpg_positions)
    if file_name != "calls":
        bias_pos = [
            (
                "normal"
                if key in tool_dict and key in true_dict
                else "not in true" if key in tool_dict else "not in tool"
            )
            for key in all_cpg_positions
        ]
        # prob_pos = [1 for _ in all_cpg_positions]
    else:
        bias_pos = [
            (
                bias_dict[key]
                if key in tool_dict and key in true_dict
                else "not in true" if key in tool_dict else "not in tool"
            )
            for key in all_cpg_positions
        ]

        # prob_pos = [
        #     prob_dict[key] if key in prob_dict else 1 for key in all_cpg_positions
        # ]

    tool_meth_values = [tool_dict.get(key, 0) for key in all_cpg_positions]
    true_meth_values = [true_dict.get(key, 0) for key in all_cpg_positions]
    # Plot TrueMeth vs toolMethod
    # Anzahl der Elemente, die größer als 0.1 sind

    output = snakemake.output["plot"]
    dist_output = snakemake.output["tool_dist"]

    # If we want to have multiple plots for different coverages we have to chose the correct one
    if isinstance(output, list):
        output = output[cov_bin]
        dist_output = dist_output[cov_bin]

    if file_name == "calls":
        file_name = "Varlociraptor"

    chromosome = [chrom for (chrom, _) in all_cpg_positions]
    position = [position for (_, position) in all_cpg_positions]
    data = pd.DataFrame(
        {
            file_name: tool_meth_values,
            "TrueMeth": true_meth_values,
            "Bias": bias_pos,
            "chromosome": chromosome,
            "position": position,
        }
    )

    filtered_data = data[data["Bias"] == "normal"]

    precision, recall, coverages = compute_precision_recall(
        file_name, cov_bin, filtered_data
    )
    save_precision_recall(
        file_name,
        coverages,
        len(filtered_data),
        precision,
        recall,
        snakemake.output["precall"],
    )

    rmse = plot_distances(
        file_name, "TrueMeth", tool_meth_values, true_meth_values, bias_pos, dist_output
    )
    plot_meth_vals(
        file_name,
        "TrueMeth",
        filtered_data,
        output,
        rmse,
    )
