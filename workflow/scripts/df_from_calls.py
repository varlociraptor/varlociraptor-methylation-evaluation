import os
import re
import pandas as pd
import sys


# Redirect standard error to snakemake log file
sys.stderr = open(snakemake.log[0], "w")


def cb_to_number(s):
    """Extracts all integers not preceded by 'E' or 'e' and returns their sum."""
    matches = re.findall(r"\d+(?=[^Ee]|$)", s)
    return sum(map(int, matches))


def get_bin(coverage):
    """Assigns a coverage bin to the given coverage value."""
    return min(
        int(coverage / snakemake.params["cov_bin_size"]),
        snakemake.params["cov_bins"] - 1,
    )


def compute_bias(format_values):
    """Returns the bias label if present in FORMAT values, else 'normal'."""
    bias_labels = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]
    if any(value != "." for value in format_values[6:12]):
        return bias_labels[format_values[6:12].index(".")]
    return "normal"


def read_simulated_file(filepath):
    """Parses simulated truth and coverage file."""
    cov_records, truth_records = [], []
    with open(filepath, "r") as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            parts1 = lines[i].strip().split("\t")
            chrom1, pos = parts1[0], (float(parts1[1]) + float(parts1[2])) / 2
            meth_reads1, cov1 = map(float, parts1[3:5])

            # if type1 == "CG" and i + 1 < len(lines):
            #     parts2 = lines[i + 1].strip().split("\t")
            #     type2 = parts2[3]
            #     if type2 == "CG":
            #         pos2 = int(parts2[1])
            #         meth_reads2, cov2 = map(float, parts2[4:6])
            #         coverage = (cov1 + cov2) / 2
            #         meth_rate = (
            #             ((meth_reads1 + meth_reads2) / (cov1 + cov2) * 100)
            #             if (cov1 + cov2)
            #             else 0
            #         )
            cov_records.append([chrom1, pos, cov1, get_bin(cov1)])
            truth_records.append([chrom1, pos, meth_reads1])
            # i += 2
            # continue
            i += 1
    return (
        pd.DataFrame(
            cov_records, columns=["chromosome", "position", "coverage", "cov_bin"]
        ),
        pd.DataFrame(
            truth_records, columns=["chromosome", "position", "true_methylation"]
        ),
    )


def read_coverage_file(filepath):
    """Reads general coverage data into a DataFrame."""
    records = []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom = parts[0]
            position = int(parts[1])
            coverage = float(parts[3])
            records.append([chrom, position, coverage, get_bin(coverage)])
    return pd.DataFrame(
        records, columns=["chromosome", "position", "coverage", "cov_bin"]
    )


def read_truth_file(filepath):
    """Reads true methylation data from BED-like format."""
    records = []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            chrom = parts[0].replace("chr", "")
            start, end = map(int, parts[1:3])
            meth_rate = float(parts[3])
            position = (start + end) // 2
            records.append([chrom, position, meth_rate])
    return pd.DataFrame(records, columns=["chromosome", "position", "true_methylation"])


def read_tool_file(filepath, file_name):
    """Reads tool-specific output files and formats them."""
    records = []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")

            if file_name == "varlo":
                chrom = parts[0]
                position = int(parts[1])
                alternative = parts[4]
                info_field = parts[7]
                format_field = parts[8]
                values = parts[9].split(":")

                if alternative != "<METH>":
                    continue

                format_fields = format_field.split(":")
                dp_index = format_fields.index("DP")
                af_index = format_fields.index("AF")
                coverage = int(values[dp_index])

                if snakemake.params["meth_type"] == "ratio":
                    saobs = cb_to_number(values[format_fields.index("SAOBS")])
                    srobs = cb_to_number(values[format_fields.index("SROBS")])
                    meth_rate = (
                        (saobs / (saobs + srobs) * 100) if (saobs + srobs) else 0
                    )
                else:
                    meth_rate = float(values[af_index]) * 100

                info_dict = dict(
                    item.split("=", 1)
                    for item in info_field.strip().split(";")
                    if "=" in item
                )

                # Zugriff auf Werte
                prob_high = info_dict["PROB_HIGH"]
                prob_low = info_dict["PROB_LOW"]
                prob_absent = info_dict["PROB_ABSENT"]


                try:
                    # TODO: Shoud I inclue PROB_LOW?
                    prob_high = 10 ** (-float(prob_high) / 10)
                    prob_low = 10 ** (-float(prob_low) / 10)
                    prob_present = prob_high + prob_low
                    prob_absent = 10 ** (-float(prob_absent) / 10)
                    if (
                        prob_present < snakemake.params["prob_pres_threshhold"]
                        and prob_absent < snakemake.params["prob_absent_threshhold"]
                    ):
                        continue
                except Exception:
                    print(
                        f"Prob present not found on chrom {chrom}, position {position}",
                        file=sys.stderr,
                    )
                    continue

                bias = compute_bias(values)
                records.append(
                    [
                        chrom,
                        position,
                        meth_rate,
                        coverage,
                        get_bin(coverage),
                        bias,
                        prob_present,
                    ]
                )

            elif file_name in {"methylDackel", "bismark"}:
                chrom = parts[0]
                start, end = int(parts[1]), int(parts[2])
                meth_rate = float(parts[3])
                coverage = int(parts[4]) + int(parts[5])
                position = (start + end) // 2
                records.append(
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

            elif file_name == "bsMap":
                if parts[0].startswith("chr"):
                    continue
                chrom = parts[0]
                position = int(parts[1])
                meth_rate = float(parts[4]) * 100
                coverage = float(parts[5])
                records.append(
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
                chrom = parts[0]
                position = int(parts[2])
                meth_rate = float(parts[3])
                coverage = int(parts[4])
                records.append(
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
                chrom = parts[0].removeprefix("chr")
                position = int(parts[2])
                details = parts[9].split()
                coverage = int(details[0])
                meth_rate = float(details[1])
                records.append(
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
                chrom = parts[0]
                position = int(parts[2])
                meth_rate = float(parts[3])
                coverage = int(parts[5])
                records.append(
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

    columns = [
        "chromosome",
        "position",
        "tool_methylation",
        "tool_coverage",
        "tool_cov_bin",
        "bias",
        "prob_present",
    ]
    return pd.DataFrame(records, columns=columns)


pd.set_option("display.max_columns", None)

tool_file = snakemake.input["tool"]
file_name = os.path.splitext(os.path.basename(tool_file))[0]

tool_df = read_tool_file(tool_file, file_name)

if snakemake.params["simulated"]:
    cov_df, truth_df = read_simulated_file(snakemake.input["true_meth"][0])
else:
    truth_df = read_truth_file(snakemake.input["true_meth"][0])
    cov_df = read_coverage_file(snakemake.input["coverage"])

df = pd.merge(truth_df, tool_df, on=["chromosome", "position"], how="inner")
df = pd.merge(df, cov_df, on=["chromosome", "position"], how="inner")
df.fillna(0, inplace=True)
df = df[df["coverage"] > 0]
df["bias"] = df["bias"].astype(str)

df.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
