import os
import re
import pandas as pd
import altair as alt
import sys
from pathlib import Path
import numpy as np

sys.stderr = open(snakemake.log[0], "w")


def point_plot(df, x, y, color, shape, x_title, y_title, height=140):

    meth_caller_to_name = {
        "varlociraptor": "Varlociraptor",
        "bismark": "Bismark",
        "bsmapz": "BSMAPz",
        "methylDackel": "MethylDackel",
        "modkit": "Modkit",
        "pb-CpG-tools": "pb-CpG-tools",
    }

    tool_base_colors = {
        "Bismark": "#D81B60",
        "BSMAPz": "#1E88E5",
        "MethylDackel": "#FFC107",
        "Varlociraptor": "#004D40",
        "Modkit": "#B42CEA",
        "pb-CpG-tools": "#8D9279",
    }

    # Map nicer names
    df = df.copy()
    df["tool_label"] = df[color].map(meth_caller_to_name)

    charts = []
    for platform in df["platform"].unique():

        subset = df[df["platform"] == platform]
        # tools that appear in this subplot
        present_tools = subset["tool_label"].unique().tolist()
        color_domain = present_tools
        color_range = [tool_base_colors[t] for t in present_tools]
        ticks = list(
            np.logspace(0, np.log10(subset[x].max() * 1.3), num=5).round().astype(int)
        )
        base = (
            alt.Chart(subset)
            .encode(
                x=alt.X(
                    f"{x}:Q",
                    title=x_title,
                    scale=alt.Scale(type="log", domain=[1, subset[x].max() * 1.3 + 1]),
                    axis=alt.Axis(values=ticks, labelAngle=-45),
                ),
                y=alt.Y(f"{y}:Q", title=y_title, scale=alt.Scale(type="log")),
                color=alt.Color(
                    "tool_label:N",
                    title="Caller",
                    scale=alt.Scale(domain=color_domain, range=color_range),
                ),
                shape=alt.Shape(f"{shape}:N") if shape else alt.value("circle"),
                tooltip=(
                    [f"{x}:Q", f"{y}:Q", "tool_label:N", f"{shape}:N"]
                    if shape
                    else [f"{x}:Q", f"{y}:Q", "tool_label:N"]
                ),
            )
            .mark_point(filled=False, size=30)
            .properties(height=height, width=height, title=platform)
            .interactive()
        )

        charts.append(base)

    return alt.hconcat(*charts).resolve_scale(color="independent")


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

df_all["task"] = np.where(
    df_all["meth_caller"] != "varlociraptor",
    "calling",
    df_all["task"],
)
# Compare different methylation calling callers
df_compare_callers = (
    df_all.groupby(["platform", "meth_caller", "replicate", "task"], as_index=False)
    .agg(s=("s", "sum"), max_rss=("max_rss", "max"))
    .query("platform != 'multi_sample'")
    .assign(
        caller=lambda x: x["meth_caller"].replace(
            {"calling": "Varlociraptor", "preprocessing": "Varlociraptor"}
        ),
        minutes=lambda x: x["s"] / 60,
        max_rss_gb=lambda x: x["max_rss"] / 1024,
    )
    .groupby(["platform", "caller", "replicate", "task"], as_index=False)
    .agg({"minutes": "sum", "max_rss_gb": "max"})
    .assign(platform_caller=lambda x: x["platform"] + " - " + x["caller"])
)


# Create runtime and memory plots for all callers
runtime_chart = point_plot(
    df_compare_callers,
    x="minutes",
    y="max_rss_gb",
    color="caller",
    shape="task",
    x_title="Runtime (min)",
    y_title="Max RSS (GB)",
    height=140,
)

# Combine caller plots
# caller_chart = alt.hconcat(runtime_chart, memory_chart)
runtime_chart.save(
    snakemake.output["tools"],
    embed_options={"actions": False},
    inline=False,
)
