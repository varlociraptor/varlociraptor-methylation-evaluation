import pandas as pd
import numpy as np
import altair as alt
import csv
import os


def compute_precision_recall(
    df, cov_bin, methylation_threshold, filename="tool", prob_threshold=0, last=False
):
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

    if prob_threshold == 0.99:
        filtered_rows = df[(df["tool_binary"] == 1) & (df["truth_binary"] == 0)]
        print("Debug")
        print(filtered_rows.to_string())
        print(FP)

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    cov_bin_size = int(snakemake.params["cov_bin_size"])

    if last:
        coverages = f"{cov_bin* cov_bin_size}-max_cov"
    elif cov_bin == "all":
        coverages = f"0 - max_cov"
    else:
        cov_bin = int(cov_bin)
        coverages = f"{cov_bin* cov_bin_size}-{cov_bin* cov_bin_size+ cov_bin_size-1}"

    return precision, recall, coverages, TP, FP, FN


def compute_precision_recall_thresholds(df, methylation_threshold, last=False):
    df = df.sort_values(by="prob_present", ascending=False).reset_index(drop=True)

    (
        precision_list,
        recall_list,
        coverages_list,
        prob_present_threshhold,
        number_sites,
        TP_list,
        FP_list,
        FN_list,
    ) = ([], [], [], [], [], [], [], [])
    # Schritt 2: Teile die Zeilen in 100 gleich große Batches (so gut wie möglich)

    # Schritt 3: Für jeden Batch etwas berechnen
    thresholds = np.arange(0, 1, 0.01)
    for prob_threshold in thresholds:
        # filtered_df = df[df["prob_present"] > threshold]
        precision, recall, coverages, TP, FP, FN = compute_precision_recall(
            df, cov_bin, methylation_threshold, "varlo", prob_threshold, last
        )
        precision_list.append(precision)
        recall_list.append(recall)
        coverages_list.append(coverages)
        prob_present_threshhold.append(prob_threshold)
        number_sites.append(len(df))
        TP_list.append(TP)
        FP_list.append(FP)
        FN_list.append(FN)
    return (
        precision_list,
        recall_list,
        coverages_list,
        prob_present_threshhold,
        number_sites,
        TP_list,
        FP_list,
        FN_list,
    )


def save_precision_recall(
    tool,
    coverage,
    no_sites,
    precision,
    recall,
    prob_present_threshhold,
    TP,
    FP,
    FN,
    filename,
):
    # Datei anlegen mit Header, wenn sie nicht existiert
    try:
        with open(filename, "x", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(
                [
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
            )
    except FileExistsError:
        pass  # Datei existiert bereits, also Header nicht neu schreiben

    # Ergebnisse anhängen
    with open(filename, "a", newline="") as file:
        writer = csv.writer(file)
        for i in range(0, len(precision)):
            writer.writerow(
                [
                    tool,
                    coverage[i],
                    no_sites[i],
                    precision[i],
                    prob_present_threshhold[i],
                    recall[i],
                    TP[i],
                    FP[i],
                    FN[i],
                ]
            )


# Einstellungen für DataFrames
alt.data_transformers.enable("vegafusion")
pd.set_option("display.max_columns", None)

# Konstanten für Berechnungen
methylation_threshold = 5  # Methylation threshold


tool_file = snakemake.input[0]
file_name = os.path.splitext(os.path.basename(tool_file))[0]

df = pd.read_parquet(tool_file, engine="pyarrow")
print("cols:", df.columns)
print(df.to_string())
max_cov = df["coverage"].max()
cov_bin = snakemake.params["cov_bin"]

if cov_bin != "all":
    df = df[df["cov_bin"] == int(cov_bin)]
print(df)

if file_name == "varlo" and snakemake.params["meth_type"] == "posterior":
    precision, recall, coverages, prob_present_threshhold, number_sites, TP, FP, FN = (
        compute_precision_recall_thresholds(
            df,
            methylation_threshold,
            # last=cov_bin == len(snakemake.output["plot"]) - 1,
        )
    )
    save_precision_recall(
        file_name,
        coverages,
        number_sites,
        precision,
        recall,
        prob_present_threshhold,
        TP,
        FP,
        FN,
        snakemake.output["precall"],
    )
else:
    precision, recall, coverages, TP, FP, FN = compute_precision_recall(
        df,
        cov_bin,
        methylation_threshold,
        # last=cov_bin == len(snakemake.output["plot"]) - 1,
    )

    save_precision_recall(
        file_name,
        [coverages],
        [len(df)],
        [precision],
        [recall],
        [0],
        [TP],
        [FP],
        [FN],
        snakemake.output["precall"],
    )
