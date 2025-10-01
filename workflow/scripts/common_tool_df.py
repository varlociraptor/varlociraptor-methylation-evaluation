import pandas as pd
import os
from functools import reduce
import altair as alt
import numpy as np
import sys


sys.stderr = open(snakemake.log[0], "w")

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)




# ---- Daten laden und zusammenführen ----
tool_dfs = []
tool_names = []
# tool_files = snakemake.input["tools"] + [snakemake.input["varlo"]]
tool_files = [snakemake.input["rep1"], snakemake.input["rep2"]]

for tool_file in tool_files:
    print(tool_file, file=sys.stderr)
    tool_name = os.path.splitext(os.path.basename(tool_file))[0]
    tool_names.append(tool_name)

    df = pd.read_parquet(tool_file, engine="pyarrow")
    print(df.head(), file=sys.stderr)
    # if tool_name == "varlo":
    #     df = df[df["bias"] == "normal"]

    # Nur notwendige Spalten behalten und umbenennen
    df = df[
        [
            "chromosome",
            "position",
            "tool_methylation",
            # "tool_coverage",
            # "prob_present",
            # "true_methylation",
        ]
    ]
    df = df.rename(
        columns={
            "tool_methylation": f"{tool_name}_methylation",
            # "tool_coverage": f"{tool_name}_coverage",
            # "prob_present": f"{tool_name}_prob_present",
        }
    )
    tool_dfs.append(df)

# Daten zusammenführen
df_merged = reduce(
    lambda left, right: pd.merge(
        left, right, on=["chromosome", "position"], how="outer"
    ),
    tool_dfs,
)
print(df_merged, file=sys.stderr)

# Ergebnis abspeichern
df_merged.to_parquet(
    snakemake.output["protocol_df"], engine="pyarrow", compression="snappy"
)
