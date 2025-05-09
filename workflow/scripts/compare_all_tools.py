import pandas as pd
import os
from functools import reduce
import altair as alt

sys.stderr = open(snakemake.log[0], "w")


def distance_plot(melted_df):
    chart = (
        alt.Chart(melted_df)
        .mark_line()
        .encode(
            x=alt.X("Distance:Q", title="Distance (%)"),
            y=alt.Y("count():Q", title="Count"),
            color=alt.Color("Tool:N", legend=alt.Legend(title="Tool")),
        )
        .properties(
            title=alt.TitleParams(
                text="Distance Comparison",
                subtitle="Comparison of the distances between predicted and true methylation rates for each tool.",
            )
        )
    )

    return chart


def compute_distances(df, tool_names):
    for tool in tool_names:
        df_short = df.dropna(subset=[f"{tool}_methylation", "true_methylation"])
        df[f"{tool}_distance"] = (
            (df_short[f"{tool}_methylation"] - df_short["true_methylation"])
            .abs()
            .round()
            .astype(int)
        )
    dist_cols = [f"{tool}_distance" for tool in tool_names]
    melted_data = df.melt(
        value_vars=dist_cols,
        var_name="Tool",
        value_name="Distance",
    )

    melted_data["Tool"] = melted_data["Tool"].str.replace("_distance", "")

    return melted_data


def density_plot(df, tool_names, target):
    melt_cols = [f"{t}_{target}" for t in tool_names]
    melted = df.melt(value_vars=melt_cols, var_name="Tool", value_name=target)
    melted["Tool"] = melted["Tool"].str.replace(f"_{target}", "")
    chart = (
        alt.Chart(melted)
        .transform_density(
            target,
            as_=[target, "Density"],
            groupby=["Tool"],
            extent=[0, 100],  # Optional: beschränkt auch die Berechnung
        )
        .transform_filter(f"datum['{target}'] <= 100")
        .mark_line()
        .encode(
            x=alt.X(
                f"{target}:Q", title=f"{target} Level", scale=alt.Scale(domain=[0, 100])
            ),
            y=alt.Y("Density:Q", title="Density"),
            color=alt.Color("Tool:N", title="Tool"),
        )
        .properties(title=f"Distribution of {target} Levels per Tool")
    )
    return chart


pd.set_option("display.max_columns", None)

tool_dfs = []
tool_names = []
for tool_file in snakemake.input["tools"]:
    tool_name = os.path.splitext(os.path.basename(tool_file))[0]
    tool_names.append(tool_name)

    df = pd.read_parquet(tool_file, engine="pyarrow")

    if tool_name == "varlo":
        df = df[df["prob_present"] >= float(snakemake.params["prob_pres_threshhold"])]
        df = df[df["bias"] == "normal"]
    df = df[
        [
            "chromosome",
            "position",
            "tool_methylation",
            "tool_coverage",
            "true_methylation",
        ]
    ].copy()
    df = df.rename(
        columns={
            "tool_methylation": f"{tool_name}_methylation",
            "tool_coverage": f"{tool_name}_coverage",
        }
    )
    tool_dfs.append(df)


df_merged = reduce(
    lambda left, right: pd.merge(
        left, right, on=["chromosome", "position", "true_methylation"], how="outer"
    ),
    tool_dfs,
)
df_merged.to_parquet(
    snakemake.output["protocol_df"], engine="pyarrow", compression="snappy"
)

methylation_density_chart = density_plot(df_merged, tool_names, "methylation")
coverage_density_chart = density_plot(df_merged, tool_names, "coverage")

melted_data = compute_distances(df_merged, tool_names)
distance_chart = distance_plot(melted_data)

chart = methylation_density_chart | coverage_density_chart | distance_chart


chart.save(snakemake.output["plot"], scale_factor=2.0)
