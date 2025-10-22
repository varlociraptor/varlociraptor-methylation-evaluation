import sys
import os
from functools import reduce
import pandas as pd
import altair as alt
import numpy as np

# Redirect stderr to Snakemake log
sys.stderr = open(snakemake.log[0], "w")

tool_dfs = []
tool_names = []

# Combine standard tool files and Varlociraptor output
tool_files = snakemake.input["tools"] + [snakemake.input["varlo"]]

for tool_file in tool_files:
    tool_name = os.path.splitext(os.path.basename(tool_file))[0]
    tool_names.append(tool_name)

    df = pd.read_parquet(tool_file, engine="pyarrow")

    # Keep only relevant columns and rename methylation column
    df = df[["chromosome", "position", "tool_methylation"]].rename(
        columns={"tool_methylation": f"{tool_name}_methylation"}
    )

    tool_dfs.append(df)

# Merge all tool data on chromosome and position
df_merged = reduce(
    lambda left, right: pd.merge(
        left, right, on=["chromosome", "position"], how="outer"
    ),
    tool_dfs,
)

df_merged.to_parquet(
    snakemake.output["protocol_df"], engine="pyarrow", compression="snappy"
)
