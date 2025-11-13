import os
import re
import pandas as pd
import altair as alt
import sys
from pathlib import Path


sys.stderr = open(snakemake.log[0], "w")


def point_plot(df, x, y, color, x_title, y_title, title, height=200):
    """Simple reusable Altair boxplot with overlaid points."""
    colorblind_safe_palette = [
        "#D81B60",
        "#1E88E5",
        "#FFC107",
        "#D35892",
        "#AC3FE6",
    ]
    base = (
        alt.Chart(df)
        .encode(
            x=alt.X(
                f"{x}:N",
                title=None,
                axis=alt.Axis(labelExpr='split(datum.value, " - ")[1]'),
            ),
            y=alt.Y(f"{y}:Q", title=y_title, scale=alt.Scale(type="log")),
            color=alt.Color(
                f"{color}:N", scale=alt.Scale(range=colorblind_safe_palette)
            ),
            xOffset=alt.XOffset(
                "jitter:Q",
                scale=alt.Scale(domain=[-0.05, 0.05]),  # << feste min/max Jitter Range
            ),
        )
        .transform_calculate(
            # jitter sehr klein!
            jitter="(random() - 0.5) * 0.02"
        )
    )

    points = base.mark_point(filled=True, size=20).encode(tooltip=[x, y])
    return points.properties(height=height)


# Read benchmark files from Snakemake input directory
records = []
benchmark_path = snakemake.input.benchmarks


for root, _, files in os.walk(benchmark_path):
    for fname in files:
        if not fname.endswith(".benchmark.txt"):
            continue

        full_path = os.path.join(root, fname)
        df = pd.read_csv(full_path, sep="\t", usecols=["s", "max_rss"])
        p = Path(full_path)

        # Extract info from folder structure
        df["platform"] = p.parts[-4]  # → "Illumina_pe"
        df["meth_caller"] = p.parts[-3]  # → "varlociraptor"
        df["task"] = p.parts[-2]  # → "simulated_data_1"
        replicate = re.sub(r"_\d+-of-\d+", "", fname.replace(".benchmark.txt", ""))
        df["replicate"] = replicate
        # Clean up sample/replicate name
        records.append(df)

df_all = pd.concat(records, ignore_index=True)
df_all["platform"] = df_all["platform"].replace("Illumina_pe", "Illumina")


# Compare different methylation calling tools
df_compare_tools = (
    df_all.groupby(["platform", "meth_caller", "replicate"], as_index=False)
    .agg(s=("s", "sum"), max_rss=("max_rss", "max"))
    .query("platform != 'multi_sample'")
    .assign(
        tool=lambda x: x["meth_caller"].replace(
            {"calling": "Varlociraptor", "preprocessing": "Varlociraptor"}
        ),
        minutes=lambda x: x["s"] / 60,
        max_rss_gb=lambda x: x["max_rss"] / 1024,
    )
    .groupby(["platform", "tool", "replicate"], as_index=False)
    .agg({"minutes": "sum", "max_rss_gb": "max"})
    .assign(platform_tool=lambda x: x["platform"] + " - " + x["tool"])
)


# Compare Varlociraptor steps individually
df_compare_varlo = (
    df_all.groupby(["platform", "meth_caller", "replicate", "task"], as_index=False)
    .agg(s=("s", "sum"), max_rss=("max_rss", "max"))
    .query("task in ['calling', 'preprocessing']")
    .assign(
        minutes=lambda x: x["s"] / 60,
        max_rss_gb=lambda x: x["max_rss"] / 1024,
        platform_tool=lambda x: x["platform"] + " - " + x["task"],
    )
)
# Create runtime and memory plots for all tools
runtime_chart = point_plot(
    df_compare_tools,
    x="platform_tool",
    y="minutes",
    color="platform",
    x_title="Tool",
    y_title="Runtime (min)",
    title="Runtime",
    height=200,
)

memory_chart = point_plot(
    df_compare_tools,
    x="platform_tool",
    y="max_rss_gb",
    color="platform",
    x_title="Platform - Tool",
    y_title="Max RSS (GB)",
    title="Memory usage",
    height=200,
)

# Combine tool plots
tool_chart = alt.hconcat(runtime_chart, memory_chart)
tool_chart.save(
    snakemake.output["tools"],
    embed_options={"actions": False},
    inline=False,
)

df_compare_varlo["platform"] = " - " + df_compare_varlo["platform"]
# Runtime comparison for Varlociraptor steps
runtime_varlo_chart = point_plot(
    df_compare_varlo,
    x="platform",
    y="minutes",
    color="task",
    x_title="Platform - Step",
    y_title="Runtime (min)",
    title="Runtime of Varlociraptor steps",
    height=150,
)

runtime_varlo_chart.save(
    snakemake.output["varlo"],
    embed_options={"actions": False},
    inline=False,
)
