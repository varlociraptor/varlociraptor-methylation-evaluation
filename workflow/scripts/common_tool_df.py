import gc
import os
import sys

import pandas as pd

# Redirect stderr to Snakemake log
sys.stderr = open(snakemake.log[0], "w")

# Combine standard tool files and Varlociraptor output
tool_files = snakemake.input["tools"] + snakemake.input["varlo"]
filter_chrom = snakemake.params["filter_chrom"]
dfs = []

for tool_file in tool_files:
    tool_name = os.path.splitext(os.path.basename(tool_file))[0]

    df = pd.read_parquet(
        tool_file,
        engine="pyarrow",
        columns=["chromosome", "position", "tool_methylation", "format"],
    )

    if filter_chrom is not None:
        df = df[df["chromosome"] == filter_chrom]
    if tool_name == "varlo":
        # Rename Varlociraptor columns to match other tools
        fdr = os.path.basename(
            os.path.dirname(os.path.dirname(os.path.dirname(tool_file)))
        )
        tool_name = f"varlo_{fdr}"

    df = df.rename(
        columns={
            "tool_methylation": f"{tool_name}_methylation",
            "format": f"{tool_name}_format",
        }
    )
    # Set chromosome and position as index before appending
    df = df.set_index(["chromosome", "position"])
    dfs.append(df)

df_merged = pd.concat(dfs, axis=1, join="outer").reset_index()
df_merged.to_parquet(
    snakemake.output["sample_df"], engine="pyarrow", compression="snappy"
)
