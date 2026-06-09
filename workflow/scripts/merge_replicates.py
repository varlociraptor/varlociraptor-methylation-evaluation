import re
import sys
from pathlib import Path

import pandas as pd

# Redirect stderr to Snakemake log file
sys.stderr = open(snakemake.log[0], "w")

# Pandas display options (useful for debugging)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)


def parse_sample_and_rep(replicate_name: str):
    """
    Extract sample name and replicate number.
    Example:
        EMSeq_HG002_LAB02_REP02 -> (EMSeq_HG002_LAB02, 2)
    """
    m = re.match(r"(.+)_REP0*(\d+)$", replicate_name)
    if m:
        return m.group(1), int(m.group(2))
    return replicate_name, None


replicate_dfs = {}

# ---- Load and organize replicates ---- #

for sample_file in snakemake.input:
    replicate_name = Path(sample_file).stem.removeprefix("sample_df_")
    df = pd.read_parquet(sample_file, engine="pyarrow")

    sample_name, rep = parse_sample_and_rep(replicate_name)

    if sample_name not in replicate_dfs:
        replicate_dfs[sample_name] = {}

    replicate_dfs[sample_name][rep] = df


merged_samples = {}

# ---- Merge replicates (inner join on genomic positions) ---- #

for sample_name, reps in replicate_dfs.items():
    if 1 not in reps or 2 not in reps:
        raise ValueError(f"Missing REP1 or REP2 for sample {sample_name}")

    df1 = reps[1]
    df2 = reps[2]

    merged_samples[sample_name] = pd.merge(
        df1,
        df2,
        on=["chromosome", "position"],
        how="inner",
        suffixes=("_rep1", "_rep2"),
    )

# ---- Combine all samples ---- #

combined_df = pd.concat(
    [df.assign(sample=sample_name) for sample_name, df in merged_samples.items()],
    ignore_index=True,
)

combined_df.to_parquet(
    snakemake.output[0],
    engine="pyarrow",
    index=False,
)
