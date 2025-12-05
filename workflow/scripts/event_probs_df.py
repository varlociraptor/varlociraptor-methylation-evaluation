import sys
import os
import re
import pandas as pd

# Redirect stderr to Snakemake log
sys.stderr = open(snakemake.log[0], "w")

# Display options for debugging
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)


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
            chrom = parts[0]
            position = int(parts[1])
            info_field = parts[7]

            # Parse INFO fields
            info_dict = dict(
                item.split("=", 1)
                for item in info_field.strip().split(";")
                if "=" in item
            )

            alpha = float(snakemake.params["alpha"])

            def phred_to_prob(score):
                return pd.NA if score == "." else 10 ** (-float(score) / 10)

            prob_present = phred_to_prob(info_dict["PROB_PRESENT"])
            prob_absent = phred_to_prob(info_dict.get("PROB_ABSENT", 0))
            prob_artifact = phred_to_prob(info_dict.get("PROB_ARTIFACT", 0))
            is_na = (
                pd.isna(prob_absent) or pd.isna(prob_artifact) or pd.isna(prob_present)
            )
            if alpha == 1.0:
                records.append(
                    [chrom, position, prob_present, prob_absent, prob_artifact]
                )
                continue

            if is_na or max(prob_present, prob_absent + prob_artifact) < (1 - alpha):
                print(
                    f"Low confidence site skipped: {chrom}:{position}",
                    file=sys.stderr,
                )
                continue
            records.append([chrom, position, prob_present, prob_absent, prob_artifact])

    return pd.DataFrame(
        records,
        columns=[
            "chromosome",
            "position",
            "prob_present",
            "prob_absent",
            "prob_artifact",
        ],
    )


# -----------------------------
# Main execution
# -----------------------------
tool_file = snakemake.input["tool"]
file_name = os.path.splitext(os.path.basename(tool_file))[0]

df = read_tool_file(tool_file, file_name)

# Save standardized output
df.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
