import pandas as pd
import os
import re


def number_reads(s):
    numbers = re.findall(r"\d+", s)
    if not numbers:
        return 0
    return sum(int(num) for num in numbers)


def get_bin(coverage):
    return min(
        int(coverage / snakemake.params["cov_bin_size"]),
        snakemake.params["cov_bins"] - 1,
    )


def compute_bias(prob, format_values):
    # min_value = min(prob_values)
    # min_index = prob_values.index(min_value)

    bias_labels = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]

    if any(value != "." for value in format_values[6:12]):
        bias = bias_labels[format_values[6:12].index(".")]
    else:
        bias = "normal"

    # if min_index != 0:
    #     bias = "normal"
    return bias


def read_truth_file(truth_file_path):
    # Lies TrueMeth-Daten
    df = []

    with open(truth_file_path, "r") as truth_file:
        for line in truth_file:
            if line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            chrom, start, end, meth_rate = (
                parts[0].replace("chr", ""),
                int(parts[1]),
                int(parts[2]),
                float(parts[3]),
            )
            position = (start + end) // 2
            df.append([chrom, position, meth_rate])

    return pd.DataFrame(df, columns=["chromosome", "position", "true_methylation"])


def read_tool_file(tool_file_path, file_name):
    # Allgemeine Struktur für die Datei
    df = []
    max_cov = 0
    with open(tool_file_path, "r") as tool_file:

        for line in tool_file:
            if (
                line.startswith("#")
                or line.startswith("track")
                or line.startswith("chr")
            ):
                continue
            parts = line.strip().split("\t")

            # Spezifische Struktur für jedes Tool
            if file_name == "varlo":  # Varlociraptor
                chrom, position, alternative, info_field, format_field, values = (
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
                    # meth_rate = float(values[af_index]) * 100
                    saobs_index = format_fields.index("SAOBS")
                    srobs_index = format_fields.index("SROBS")
                    saobs = number_reads(values[saobs_index])
                    srobs = number_reads(values[srobs_index])
                    if srobs == 0:
                        meth_rate = 100
                    else:
                        meth_rate = (saobs / (saobs + srobs)) * 100
                    coverage = int(values[dp_index])

                    match = re.search(r"PROB_HIGH=([\d\.]+)", info_field)
                    prob_present = 0
                    try:
                        prob_present = 10 ** (-float(match.group(1)) / 10)
                        if position == 5030481:
                            print(prob_present)
                    except:
                        print(
                            f"Prob present not found on chrom {chrom}, position {position}"
                        )
                        continue

                    bias = compute_bias(info_field, values)
                    cov_bin = get_bin(coverage)
                    if position == 5030481:
                        print(bias)
                    df.append(
                        [
                            chrom,
                            position,
                            meth_rate,
                            coverage,
                            cov_bin,
                            bias,
                            prob_present,
                        ]
                    )

            elif file_name == "methylDackel":
                chrom, start, end, meth_rate, coverage = (
                    parts[0],
                    int(parts[1]),
                    int(parts[2]),
                    float(parts[3]),
                    int(parts[4]) + int(parts[5]),
                )
                position = (start + end) // 2
                cov_bin = get_bin(coverage)
                df.append(
                    [chrom, position, meth_rate, coverage, cov_bin, "normal", None]
                )

            elif file_name == "bsMap":
                chrom, position, meth_rate, coverage = (
                    parts[0],
                    int(parts[1]),
                    float(parts[4]) * 100,
                    float(parts[5]),
                )
                df.append(
                    [
                        chrom,
                        position,
                        meth_rate,
                        coverage,
                        get_bin(coverage),
                        "normal",
                        None,
                    ]
                )

            elif file_name == "bismark":
                chrom, start, end, meth_rate, coverage = (
                    parts[0],
                    int(parts[1]),
                    int(parts[2]),
                    float(parts[3]),
                    int(parts[4]) + int(parts[5]),
                )
                position = (start + end) // 2
                df.append(
                    [
                        chrom,
                        position,
                        meth_rate,
                        coverage,
                        get_bin(coverage),
                        "normal",
                        None,
                    ]
                )

            elif file_name == "bisSNP":
                chrom, position, meth_rate, coverage = (
                    parts[0],
                    int(parts[2]),
                    float(parts[3]),
                    int(parts[4]),
                )
                df.append(
                    [
                        chrom,
                        position,
                        meth_rate,
                        coverage,
                        get_bin(coverage),
                        "normal",
                        None,
                    ]
                )

            elif file_name == "modkit":
                details = parts[9].split()
                chrom, position, coverage, meth_rate = (
                    parts[0].removeprefix("chr"),
                    int(parts[2]),
                    int(details[0]),
                    float(details[1]),
                )

                df.append(
                    [
                        chrom,
                        position,
                        meth_rate,
                        coverage,
                        get_bin(coverage),
                        "normal",
                        None,
                    ]
                )

            elif file_name == "pb_CpG_tools":
                chrom, position, meth_rate, coverage = (
                    parts[0],
                    int(parts[2]),
                    float(parts[3]),
                    int(parts[5]),
                )
                if coverage > max_cov:
                    max_cov = coverage

                df.append(
                    [
                        chrom,
                        position,
                        meth_rate,
                        coverage,
                        get_bin(coverage),
                        "normal",
                        None,
                    ]
                )
    # Erstelle DataFrame
    columns = [
        "chromosome",
        "position",
        "tool_methylation",
        "coverage",
        "cov_bin",
        "bias",
        "prob_present",
    ]
    return pd.DataFrame(df, columns=columns)


# Globale DataFrames
tool_df = pd.DataFrame()
truth_df = pd.DataFrame()


tool_file = snakemake.input["tool"]
file_name = os.path.splitext(os.path.basename(tool_file))[0]

tool_df = read_tool_file(tool_file, file_name)
truth_df = read_truth_file(snakemake.input["true_meth"][0])

max_cov = tool_df["coverage"].max()
print("Debug")
print(tool_df.shape)
print(tool_df[tool_df["position"] == "5030481"])
df = pd.merge(
    tool_df[tool_df["bias"] == "normal"],  # Filter Tool-Daten
    truth_df,  # Filter True-Daten
    on=["chromosome", "position"],
)
df.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
