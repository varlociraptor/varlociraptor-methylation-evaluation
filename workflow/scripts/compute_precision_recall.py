import pandas as pd
import numpy as np
import altair as alt
import csv
import os


# Redirect standard error to snakemake log file
sys.stderr = open(snakemake.log[0], "w")


# Set up plotting and display
alt.data_transformers.enable("vegafusion")
pd.set_option("display.max_columns", None)

# Threshold for methylation (When do we say a site is methylated. 0 is bad since the truth (avg bedgraph) has nearly never meth rates of 0)
methylation_threshold = 5


def compute_precision_recall(
    df, cov_bin, methylation_threshold, filename="tool", prob_threshold=0, last=False
):
    """Compute precision and recall based on thresholds and coverage bins."""
    if filename == "varlo" and snakemake.params["meth_type"] == "posterior":
        df["tool_binary"] = df["prob_present"].apply(
            lambda x: 1 if x > prob_threshold else 0
        )
    else:
        df["tool_binary"] = df["tool_methylation"].apply(
            lambda x: 1 if x > methylation_threshold else 0
        )

    df["truth_binary"] = df["true_methylation"].apply(
        lambda x: 1 if x > methylation_threshold else 0
    )

    TP = ((df["tool_binary"] == 1) & (df["truth_binary"] == 1)).sum()
    FP = ((df["tool_binary"] == 1) & (df["truth_binary"] == 0)).sum()
    FN = ((df["tool_binary"] == 0) & (df["truth_binary"] == 1)).sum()

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0

    cov_bin_size = int(snakemake.params["cov_bin_size"])

    if last:
        coverage_str = f"{cov_bin * cov_bin_size}-max_cov"
    elif cov_bin == "all":
        coverage_str = "0 - max_cov"
    else:
        cov_bin = int(cov_bin)
        start = cov_bin * cov_bin_size
        end = start + cov_bin_size - 1
        coverage_str = f"{start}-{end}"

    return precision, recall, coverage_str, TP, FP, FN


def compute_precision_recall_thresholds(df, methylation_threshold, last=False):
    """Sweep through probability thresholds and compute precision/recall at each step."""
    df = df.sort_values(by="prob_present", ascending=False).reset_index(drop=True)

    results = {
        "precision": [],
        "recall": [],
        "coverage": [],
        "thresholds": [],
        "num_sites": [],
        "TP": [],
        "FP": [],
        "FN": [],
    }

    for prob_threshold in np.arange(0, 1, 0.01):
        precision, recall, coverage_str, TP, FP, FN = compute_precision_recall(
            df,
            cov_bin,
            methylation_threshold,
            filename="varlo",
            prob_threshold=prob_threshold,
            last=last,
        )
        results["precision"].append(precision)
        results["recall"].append(recall)
        results["coverage"].append(coverage_str)
        results["thresholds"].append(prob_threshold)
        results["num_sites"].append(len(df))
        results["TP"].append(TP)
        results["FP"].append(FP)
        results["FN"].append(FN)

    return (
        results["precision"],
        results["recall"],
        results["coverage"],
        results["thresholds"],
        results["num_sites"],
        results["TP"],
        results["FP"],
        results["FN"],
    )


def save_precision_recall(
    tool, coverage, num_sites, precision, recall, thresholds, TP, FP, FN, filename
):
    """Write precision/recall results to a CSV file."""
    header = [
        "tool",
        "coverage",
        "number_sites",
        "precision",
        "prob_pres_threshold",
        "recall",
        "TP",
        "FP",
        "FN",
    ]

    # Write header only if the file is new
    try:
        with open(filename, "x", newline="") as file:
            csv.writer(file).writerow(header)
    except FileExistsError:
        pass

    # Append data
    with open(filename, "a", newline="") as file:
        writer = csv.writer(file)
        for i in range(len(precision)):
            writer.writerow(
                [
                    tool,
                    coverage[i],
                    num_sites[i],
                    precision[i],
                    thresholds[i],
                    recall[i],
                    TP[i],
                    FP[i],
                    FN[i],
                ]
            )


# Load data and determine coverage bin
tool_file = snakemake.input[0]
file_name = os.path.splitext(os.path.basename(tool_file))[0]

df = pd.read_parquet(tool_file, engine="pyarrow")
cov_bin = snakemake.params["cov_bin"]

if cov_bin != "all":
    df = df[df["cov_bin"] == int(cov_bin)]

# Compute and save precision/recall data
if file_name == "varlo" and snakemake.params["meth_type"] == "posterior":
    precision, recall, coverage, thresholds, num_sites, TP, FP, FN = (
        compute_precision_recall_thresholds(df, methylation_threshold)
    )
    save_precision_recall(
        tool=file_name,
        coverage=coverage,
        num_sites=num_sites,
        precision=precision,
        recall=recall,
        thresholds=thresholds,
        TP=TP,
        FP=FP,
        FN=FN,
        filename=snakemake.output["precall"],
    )
else:
    precision, recall, coverage_str, TP, FP, FN = compute_precision_recall(
        df, cov_bin, methylation_threshold
    )
    save_precision_recall(
        tool=file_name,
        coverage=[coverage_str],
        num_sites=[len(df)],
        precision=[precision],
        recall=[recall],
        thresholds=[0],
        TP=[TP],
        FP=[FP],
        FN=[FN],
        filename=snakemake.output["precall"],
    )
