import os
import re
import pandas as pd
import altair as alt

sys.stderr = open(snakemake.log[0], "w")


def boxplot_with_points(df, x, y, color, x_title, y_title, title):
    """Simple reusable Altair boxplot with overlaid points."""
    base = alt.Chart(df).encode(
        x=alt.X(
            f"{x}:N",
            title=x_title,
            axis=alt.Axis(labelExpr='split(datum.value, " - ")[1]'),
        ),
        y=alt.Y(f"{y}:Q", title=y_title, scale=alt.Scale(type="log")),
        color=f"{color}:N",
    )
    box = base.mark_boxplot(extent="min-max")
    points = base.mark_point(filled=True, size=50).encode(tooltip=[x, y])
    return (box + points).properties(width=600, height=400, title=title)


# Read benchmark files from Snakemake input directory
records = []
benchmark_path = snakemake.input.benchmarks
from pathlib import Path


for root, _, files in os.walk(benchmark_path):
    for fname in files:
        if not fname.endswith(".benchmark.txt"):
            continue

        full_path = os.path.join(root, fname)
        df = pd.read_csv(full_path, sep="\t", usecols=["s", "max_rss"])
        p = Path(full_path)

        # Extract info from folder structure
        df["platform"] = p.parts[1]  # → "Illumina_pe"
        df["meth_caller"] = p.parts[2]  # → "varlociraptor"
        df["task"] = p.parts[3]  # → "simulated_data_1"
        replicate = re.sub(r"_\d+-of-\d+", "", fname.replace(".benchmark.txt", ""))
        df["replicate"] = replicate
        # Clean up sample/replicate name
        records.append(df)

df_all = pd.concat(records, ignore_index=True)
# Aggregate runtime and memory usage per replicate

df_summary = df_all.groupby(
    ["platform", "meth_caller", "replicate"], as_index=False
).agg(s=("s", "sum"), max_rss=("max_rss", "max"))
print(df_summary.head())
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
print(df_all.head())

df_summary = df_all.groupby(
    ["platform", "meth_caller", "replicate", "task"], as_index=False
).agg(s=("s", "sum"), max_rss=("max_rss", "max"))

print(df_summary.to_string())
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
print(df_compare_varlo["platform_tool"])

# Print summary table for debugging/logging
# Create runtime and memory plots for all tools
runtime_chart = boxplot_with_points(
    df_compare_tools,
    x="platform_tool",
    y="minutes",
    color="platform",
    x_title="Tool",
    y_title="Runtime (min)",
    title="Runtime of methylation calling tools",
)

memory_chart = boxplot_with_points(
    df_compare_tools,
    x="platform_tool",
    y="max_rss_gb",
    color="platform",
    x_title="Platform - Tool",
    y_title="Max RSS (GB)",
    title="Memory usage of methylation calling tools",
)

# Combine tool plots
tool_chart = alt.hconcat(runtime_chart, memory_chart)
tool_chart.save(
    snakemake.output["tools"],
    embed_options={"actions": False},
    inline=False,
)

# Runtime comparison for Varlociraptor steps
runtime_varlo_chart = boxplot_with_points(
    df_compare_varlo,
    x="platform_tool",
    y="minutes",
    color="platform",
    x_title="Platform - Step",
    y_title="Runtime (min)",
    title="Runtime of Varlociraptor steps",
)

runtime_varlo_chart.save(
    snakemake.output["varlo"],
    embed_options={"actions": False},
    inline=False,
)
