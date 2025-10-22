import sys
import os
import re
import pandas as pd

# Redirect stderr to Snakemake log
sys.stderr = open(snakemake.log[0], "w")

# Display options for debugging
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)


# -----------------------------
# Helper functions
# -----------------------------


def cb_to_number(s: str) -> int:
    """
    Extract all integers from a string (ignoring 'E'/'e') and return their sum.
    """
    matches = re.findall(r"\d+(?=[^Ee]|$)", s)
    return sum(map(int, matches))


def read_tool_file(filepath: str, file_name: str) -> pd.DataFrame:
    """
    Parse tool-specific methylation output into a standardized DataFrame.

    Parameters
    ----------
    filepath : str
        Path to the tool output file.
    file_name : str
        Tool identifier (e.g., 'varlo', 'methylDackel', 'bsMap').

    Returns
    -------
    pd.DataFrame
        Columns: chromosome, position, tool_methylation
    """
    records = []

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
                info_field = parts[7]
                format_field = parts[8]

                if alt != "<METH>":
                    continue

                # Parse sample AF values
                format_fields = format_field.split(":")
                try:
                    af_index = format_fields.index("AF")
                except ValueError:
                    continue

                sample_afs = []
                for sample_fmt in parts[9:]:
                    values = sample_fmt.split(":")
                    try:
                        af = float(values[af_index])
                        sample_afs.append(af * 100)
                    except (ValueError, IndexError):
                        pass

                meth_rate = sum(sample_afs) / len(sample_afs) if sample_afs else 0

                # Parse INFO fields
                info_dict = dict(
                    item.split("=", 1)
                    for item in info_field.strip().split(";")
                    if "=" in item
                )

                alpha = float(snakemake.params["alpha"])

                def phred_to_prob(score):
                    return 10 ** (-float(score) / 10)

                def is_missing(val):
                    return val is None or val == "."

                # Determine probability of methylation
                if is_missing(info_dict.get("PROB_PRESENT")):
                    if is_missing(info_dict.get("PROB_HIGH")) or is_missing(
                        info_dict.get("PROB_LOW")
                    ):
                        print(
                            f"Missing probability info at {chrom}:{position}",
                            file=sys.stderr,
                        )
                        continue
                    prob_present = phred_to_prob(
                        info_dict["PROB_HIGH"]
                    ) + phred_to_prob(info_dict["PROB_LOW"])
                else:
                    prob_present = phred_to_prob(info_dict["PROB_PRESENT"])

                prob_absent = phred_to_prob(info_dict.get("PROB_ABSENT", 0))
                prob_artifact = phred_to_prob(info_dict.get("PROB_ARTIFACT", 0))

                if max(prob_present, prob_absent + prob_artifact) < (1 - alpha):
                    print(
                        f"Low confidence site skipped: {chrom}:{position}",
                        file=sys.stderr,
                    )
                    continue

                records.append([chrom, position, meth_rate])

            # -----------------------------
            # MethylDackel / Bismark formats
            # -----------------------------
            elif file_name in {"methylDackel", "bismark"}:
                chrom = parts[0]
                start, end = int(parts[1]), int(parts[2])
                meth_rate = float(parts[3])
                position = (start + end) // 2
                records.append([chrom, position, meth_rate])

            # -----------------------------
            # bsMap format
            # -----------------------------
            elif file_name == "bsMap":
                if parts[0].startswith("chr"):
                    continue
                chrom = parts[0]
                position = int(parts[1])
                meth_rate = float(parts[4]) * 100
                records.append([chrom, position, meth_rate])

            # -----------------------------
            # bisSNP format
            # -----------------------------
            elif file_name == "bisSNP":
                chrom = parts[0]
                position = int(parts[2])
                meth_rate = float(parts[3])
                records.append([chrom, position, meth_rate])

            # -----------------------------
            # modkit format
            # -----------------------------
            elif file_name == "modkit":
                chrom = parts[0].removeprefix("chr")
                position = int(parts[2])
                modified_base = parts[3]
                details = parts[9].split()
                meth_rate = float(details[1])

                if modified_base == "m":
                    records.append([chrom, position, meth_rate])

            # -----------------------------
            # pb_CpG_tools format
            # -----------------------------
            elif file_name == "pb_CpG_tools":
                chrom = parts[0]
                position = int(parts[2])
                meth_rate = float(parts[3])
                records.append([chrom, position, meth_rate])

    return pd.DataFrame(records, columns=["chromosome", "position", "tool_methylation"])


# -----------------------------
# Main execution
# -----------------------------
tool_file = snakemake.input["tool"]
file_name = os.path.splitext(os.path.basename(tool_file))[0]

df = read_tool_file(tool_file, file_name)

# Save standardized output
df.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
