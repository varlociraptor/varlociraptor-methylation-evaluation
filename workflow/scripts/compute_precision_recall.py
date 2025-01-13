import pandas as pd
import numpy as np
import altair as alt
import csv
import os


def compute_precision_recall(df, cov_bin, methylation_threshold, last=False):

    df["tool_binary"] = df["tool_methylation"].apply(
        lambda x: 1 if x > methylation_threshold else 0
    )
    df["truth_binary"] = df["true_methylation"].apply(
        lambda x: 1 if x > methylation_threshold else 0
    )

    TP = ((df["tool_binary"] == 1) & (df["truth_binary"] == 1)).sum()
    FP = ((df["tool_binary"] == 1) & (df["truth_binary"] == 0)).sum()
    FN = ((df["tool_binary"] == 0) & (df["truth_binary"] == 1)).sum()

    print(df)
    print(TP, FP, FN)

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    cov_bin_size = int(snakemake.params["cov_bin_size"])

    if last:
        coverages = f"{cov_bin* cov_bin_size}-max_cov"
    else:
        coverages = f"{cov_bin* cov_bin_size}-{cov_bin* cov_bin_size+ cov_bin_size-1}"

    return precision, recall, coverages


def compute_precision_recall_thresholds(df, methylation_threshold, last=False):
    df = df.sort_values(by="prob_present", ascending=False).reset_index(drop=True)

    (
        precision_list,
        recall_list,
        coverages_list,
        prob_present_threshhold,
        number_sites,
    ) = (
        [],
        [],
        [],
        [],
        [],
    )
    # Schritt 2: Teile die Zeilen in 100 gleich große Batches (so gut wie möglich)

    # Schritt 3: Für jeden Batch etwas berechnen
    thresholds = np.arange(0, 1, 0.01)
    for threshold in thresholds:
        filtered_df = df[df["prob_present"] > threshold]
        precision, recall, coverages = compute_precision_recall(
            filtered_df, cov_bin, methylation_threshold, last
        )
        precision_list.append(precision)
        recall_list.append(recall)
        coverages_list.append(coverages)
        prob_present_threshhold.append(threshold)
        number_sites.append(len(filtered_df))
    return (
        precision_list,
        recall_list,
        coverages_list,
        prob_present_threshhold,
        number_sites,
    )


def save_precision_recall(
    tool, coverage, no_sites, precision, recall, prob_present_threshhold, filename
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
                ]
            )


# Einstellungen für DataFrames
alt.data_transformers.enable("vegafusion")
pd.set_option("display.max_columns", None)

# Konstanten für Berechnungen
methylation_threshold = 5  # Methylation threshold


tool_file = snakemake.input[0]
base_name = os.path.splitext(os.path.basename(tool_file))[0]
file_name = "Varlociraptor" if base_name == "calls" else base_name

df = pd.read_parquet(tool_file, engine="pyarrow")

cov_bin = int(snakemake.params["cov_bin"])
max_cov = df["coverage"].max()


df = df[df["cov_bin"] == cov_bin]
if file_name == "varlo":
    precision, recall, coverages, prob_present_threshhold, number_sites = (
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
        snakemake.output["precall"],
    )
else:
    precision, recall, coverages = compute_precision_recall(
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
        snakemake.output["precall"],
    )
