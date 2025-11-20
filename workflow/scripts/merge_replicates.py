import sys
import re
import pandas as pd

# Redirect stderr to Snakemake log file
sys.stderr = open(snakemake.log[0], "w")

# Pandas display options (useful for debugging)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)


def normalize_sample_name(replicate_name: str) -> str:
    """
    Normalize replicate names to a sample identifier.
    """
    name = re.sub(r"_REP\d+$", "", replicate_name)  # Illumina pattern
    name = re.sub(r"REP\d+$", "REP", name)  # PacBio/ Nanopore / multi-sample pattern
    name = re.sub(r"minimal\d+$", "minimal", name)  # minimal toy example
    return name or replicate_name


replicate_dfs = {}

# ---- Merge replicates for each sample ---- #

for sample_file in snakemake.input:
    # Extract replicate name from file path
    replicate_name = (
        sample_file.split("/")[-1].replace(".parquet", "").replace("sample_df_", "", 1)
    )
    df = pd.read_parquet(sample_file, engine="pyarrow")

    # Normalize sample name (combine replicates)
    sample_name = normalize_sample_name(replicate_name)
    if sample_name in replicate_dfs:
        # Merge replicate 1 and replicate 2 data
        replicate_dfs[sample_name] = pd.merge(
            replicate_dfs[sample_name],
            df,
            on=["chromosome", "position"],
            how="inner",
            suffixes=("_rep1", "_rep2"),
        )
    else:
        replicate_dfs[sample_name] = df

combined_df = pd.concat(
    [df.assign(replicate=key) for key, df in replicate_dfs.items()],
    ignore_index=True,
)
combined_df.to_parquet(
    snakemake.output[0],
    engine="pyarrow",
    index=False,
)
