import pandas as pd
import os
import re


def cb_to_number(s):
    # Finde alle Zahlen, aber ignoriere die, die vor einem "E" oder "e" stehen
    matches = re.findall(r"\d+(?=[^Ee]|$)", s)
    return sum(map(int, matches))


def get_bin(coverage):
    return min(
        int(coverage / snakemake.params["cov_bin_size"]),
        snakemake.params["cov_bins"] - 1,
    )


def compute_bias(prob, format_values):
    bias_labels = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]

    if any(value != "." for value in format_values[6:12]):
        bias = bias_labels[format_values[6:12].index(".")]
    else:
        bias = "normal"

    return bias


def read_simulated_file(simulated_file_path):
    cov_df = []
    truth_df = []
    with open(simulated_file_path, "r") as truth_file:
        lines = truth_file.readlines()
        i = 0
        while i < len(lines):
            # Zeile aufteilen
            parts1 = lines[i].strip().split("\t")
            chrom1, type1, meth_reads1, cov1, meth_rate1 = (
                parts1[0],
                parts1[3],
                float(parts1[4]),
                float(parts1[5]),
                float(parts1[6]),
            )

            if type1 == "CG":
                # Naechste Zeile
                parts2 = lines[i + 1].strip().split("\t")
                pos2, type2, meth_reads2, cov2, meth_rate2 = (
                    int(parts2[1]),
                    parts2[3],
                    float(parts2[4]),
                    float(parts2[5]),
                    float(parts2[6]),
                )

                if type2 == "CG":
                    # Aggregiere die Daten
                    position = pos2  # Mittelwert der Position
                    coverage = (cov1 + cov2) / 2  # Durchschnitt der Coverage
                    meth_rate = (
                        ((meth_reads1 + meth_reads2) / (cov1 + cov2)) * 100
                        if (cov1 + cov2) != 0
                        else 0
                    )
                    # Durchschnitt der Methylierungsrate

                    cov_df.append([chrom1, position, coverage, get_bin(coverage)])
                    truth_df.append([chrom1, position, meth_rate])
                    i += 2  # Weiter zur nächsten Zeile nach dem Paar
                else:
                    i += 1  # Falls die naechste Zeile nicht CG ist, weiter zur nächsten
            else:
                i += 1  # Falls aktuelle Zeile kein CG ist, weiter zur nächsten

    return (
        pd.DataFrame(cov_df, columns=["chromosome", "position", "coverage", "cov_bin"]),
        pd.DataFrame(truth_df, columns=["chromosome", "position", "true_methylation"]),
    )


def read_coverage_file(file_path):
    df = []

    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            chrom = parts[0]  # Chromosom
            position = int(parts[1])  # Startposition
            coverage = float(parts[3])  # Coverage-Wert
            cov_bin = get_bin(coverage)

            df.append([chrom, position, coverage, cov_bin])

    return pd.DataFrame(df, columns=["chromosome", "position", "coverage", "cov_bin"])


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
                    if snakemake.params["meth_type"] == "ratio":
                        saobs_index = format_fields.index("SAOBS")
                        srobs_index = format_fields.index("SROBS")
                        saobs = cb_to_number(values[saobs_index])
                        srobs = cb_to_number(values[srobs_index])
                        meth_rate = (
                            (saobs / (saobs + srobs) * 100)
                            if (saobs + srobs) != 0
                            else 0
                        )

                    else:
                        meth_rate = float(values[af_index]) * 100
                    coverage = int(values[dp_index])

                    match = re.search(r"PROB_HIGH=([\d\.]+)", info_field)
                    prob_present = 0
                    try:
                        prob_present = 10 ** (-float(match.group(1)) / 10)
                    except:
                        print(
                            f"Prob present not found on chrom {chrom}, position {position}"
                        )
                        if position == 5030481:
                            print("CONTINUIE")
                        continue

                    bias = compute_bias(info_field, values)
                    cov_bin = get_bin(coverage)
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
                if line.startswith("chr"):
                    continue
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
                        0,
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
                        0,
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
                        0,
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
                        0,
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
                        0,
                    ]
                )
    # Erstelle DataFrame
    columns = [
        "chromosome",
        "position",
        "tool_methylation",
        "tool_coverage",
        "tool_cov_bin",
        "bias",
        "prob_present",
    ]
    return pd.DataFrame(df, columns=columns)


pd.set_option("display.max_columns", None)
# Globale DataFrames
tool_df = pd.DataFrame()
truth_df = pd.DataFrame()


tool_file = snakemake.input["tool"]
file_name = os.path.splitext(os.path.basename(tool_file))[0]

tool_df = read_tool_file(tool_file, file_name)


if snakemake.params["simulated"]:
    cov_df, truth_df = read_simulated_file(snakemake.input["simulated"][0])
else:
    truth_df = read_truth_file(snakemake.input["true_meth"][0])
    cov_df = read_coverage_file(snakemake.input["coverage"])

# Mergen von tool_df mit truth_df (Left Join → Alle Positionen aus truth_df bleiben erhalten)
df = pd.merge(truth_df, tool_df, on=["chromosome", "position"], how="left")

# Mergen von df mit cov_df (wieder Left Join → Nur Positionen aus truth_df bleiben erhalten)
df = pd.merge(df, cov_df, on=["chromosome", "position"], how="left")

# Fehlende Werte mit 0 ersetzen
df.fillna(0, inplace=True)

print(df.columns)
print(df.shape)
print(df.head())

# Delete all values that are not covered by the reads
df = df[df["coverage"] > 0]
print(df.shape)

print(df.head())
df["bias"] = df["bias"].astype(str)
df.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
