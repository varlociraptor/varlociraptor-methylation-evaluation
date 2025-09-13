import pandas as pd
import os
from functools import reduce
import altair as alt
import numpy as np
import sys


sys.stderr = open(snakemake.log[0], "w")

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)

# stderr für Snakemake-Logging
# def distance_plot(melted_agg):
#     """Erstellt den Distance-Plot aus aggregierten Daten."""
#     chart = (
#         alt.Chart(melted_agg)
#         .mark_line()
#         .encode(
#             x=alt.X("Distance:Q", title="Distance (%)"),
#             y=alt.Y("count:Q", title="Count"),
#             color=alt.Color("Tool:N", legend=alt.Legend(title="Tool")),
#         )
#         .properties(title="Distance Comparison")
#     )
#     return chart


# def compute_distances(df, tool_names):
#     """Berechnet die gerundeten Abstände zwischen gemessener und wahrer Methylierung."""
#     for tool in tool_names:
#         df_short = df.dropna(subset=[f"{tool}_methylation", "true_methylation"])
#         df[f"{tool}_distance"] = (
#             (df_short[f"{tool}_methylation"] - df_short["true_methylation"])
#             .abs()
#             .round()
#             .astype(int)
#         )
#     dist_cols = [f"{tool}_distance" for tool in tool_names]
#     melted_data = df.melt(
#         value_vars=dist_cols,
#         var_name="Tool",
#         value_name="Distance",
#     )
#     melted_data["Tool"] = melted_data["Tool"].str.replace("_distance", "")
#     return melted_data


# # VegaFusion aktivieren, um große Datensätze effizient zu transformieren
# # alt.data_transformers.enable("vegafusion")
# def density_plot_agg(df, tool_names, target, n_bins=200):
#     """Erstellt Dichteplots durch Voraggregation in Pandas (Histogramm)."""
#     melt_cols = [f"{t}_{target}" for t in tool_names]
#     melted = df.melt(value_vars=melt_cols, var_name="Tool", value_name=target)
#     melted["Tool"] = melted["Tool"].str.replace(f"_{target}", "")

#     dfs = []
#     for tool in tool_names:
#         values = melted[melted["Tool"] == tool][target].values
#         hist, bins = np.histogram(values, bins=n_bins, range=(0, 100), density=True)
#         df_hist = pd.DataFrame({
#             target: (bins[:-1] + bins[1:]) / 2,
#             "Density": hist,
#             "Tool": tool
#         })
#         dfs.append(df_hist)
#     df_density = pd.concat(dfs, ignore_index=True)

#     chart = (
#         alt.Chart(df_density)
#         .mark_line()
#         .encode(
#             x=alt.X(f"{target}:Q", title=f"{target} Level", scale=alt.Scale(domain=[0, 100])),
#             y=alt.Y("Density:Q", title="Density"),
#             color=alt.Color("Tool:N", title="Tool")
#         )
#         .properties(title=f"Distribution of {target} Levels per Tool")
#     )
#     return chart


# ---- Daten laden und zusammenführen ----
tool_dfs = []
tool_names = []

for tool_file in snakemake.input["tools"]:
    tool_name = os.path.splitext(os.path.basename(tool_file))[0]
    tool_names.append(tool_name)

    df = pd.read_parquet(tool_file, engine="pyarrow")

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

# Ergebnis abspeichern
df_merged.to_parquet(
    snakemake.output["protocol_df"], engine="pyarrow", compression="snappy"
)
