import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)


def is_missing(val):
    return val is None or val == "."


def varlo_meth_rate_controlled(samples: list, skipped: Counter):
    """Methylation rate (%) of a site already filtered by varlociraptor's FDR control, NaN for a kept artifact.

    The calls come from `varlociraptor filter-calls control-fdr` (rule control_fdr_calls), so the site is not
    filtered again. AF as reported by Varlociraptor, averaged over the samples. Sites without AFD (artifacts kept
    by --smart-retain-artifacts) get no rate (NaN), so that they stay available for the bias plots but are not
    counted as 0 %. Returns None if no sample has an AF.
    """
    # Artifacts have no AFD and no estimable methylation rate (AF = 0 would count them as unmethylated)
    if samples and any(is_missing(values.get("AFD")) for values in samples):
        return np.nan

    sample_afs = []
    for values in samples:
        try:
            sample_afs.append(float(values["AF"]) * 100)
        except (ValueError, KeyError):
            pass
    if not sample_afs:
        skipped["missing AF"] += 1
        return None
    return sum(sample_afs) / len(sample_afs)


def read_tool_file(filepath: str, file_name: str) -> pd.DataFrame:
    """    Parse tool-specific methylation output into a standardized DataFrame.
    """
    records = []
    # Number of skipped Varlociraptor sites per reason
    skipped = Counter()

    with open(filepath, "r") as f:
        for line in f:
            # Skip header/comment lines
            if line.startswith("#") or line.startswith("track"):
                continue

            parts = line.strip().split("\t")

            # -----------------------------
            # Varlociraptor format
            # -----------------------------
            if file_name == "varlo":
                chrom = parts[0]
                position = int(parts[1])
                alt = parts[4]

                if alt != "<METH>":
                    continue

                # One FORMAT dict per sample
                format_fields = parts[8].split(":")
                samples = [
                    dict(zip(format_fields, sample_fmt.split(":")))
                    for sample_fmt in parts[9:]
                ]

                meth_rate = varlo_meth_rate_controlled(samples, skipped)
                if meth_rate is None:
                    continue
                records.append([chrom, position, meth_rate, parts[9]])

            # MethylDackel / Bismark
            elif file_name in {"methylDackel"}:
                chrom = parts[0]
                start, end = int(parts[1]), int(parts[2])
                meth_rate = float(parts[3])
                position = (start + end) // 2
                records.append([chrom, position, meth_rate, pd.NA])
            # BisSNP / Bismark
            elif file_name in {"bismark", "bisSNP"}:
                chrom = parts[0]
                position = int(parts[1])
                meth_rate = float(parts[2])
                records.append([chrom, position, meth_rate, pd.NA])

            # bsMap
            elif file_name == "bsMap":
                if parts[0].startswith("chr"):
                    continue
                chrom = parts[0]
                position = int(parts[1])
                meth_rate = float(parts[4]) * 100
                records.append([chrom, position, meth_rate, pd.NA])

            # modkit
            elif file_name == "modkit":
                chrom = parts[0].removeprefix("chr")
                position = int(parts[2])
                modified_base = parts[3]
                meth_rate = float(parts[10])

                if modified_base == "m":
                    records.append([chrom, position, meth_rate, pd.NA])

            # pb_CpG_tools
            elif file_name == "pb_CpG_tools":
                chrom = parts[0]
                position = int(parts[2])
                meth_rate = float(parts[3])
                records.append([chrom, position, meth_rate, pd.NA])

    return pd.DataFrame(
        records, columns=["chromosome", "position", "tool_methylation", "format"]
    )



tool_file = snakemake.input["tool"]
file_name = os.path.splitext(os.path.basename(tool_file))[0]

df = read_tool_file(tool_file, file_name)
df.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
