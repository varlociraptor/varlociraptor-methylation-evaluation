import os
import re
import sys
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)


def point_plot(df, x, y, color, shape, x_title, y_title, height=140):

    meth_caller_to_name = {
        "varlo_0.1": "Varlociraptor α = 0.1",
        "varlo_0.05": "Varlociraptor α = 0.05",
        "varlo_0.01": "Varlociraptor α = 0.01",
        "varlociraptor": "Varlociraptor",
        "bismark": "Bismark",
        "bsmap": "BSMAPz",
        "methylDackel": "MethylDackel",
        "modkit": "Modkit",
        "pb-CpG-tools": "Pb-CpG-tools",
        "bisSNP": "BisSNP",
    }

    tool_colors = {
        "Bismark": "#D81B60",
        "BSMAPz": "#1E88E5",
        "BisSNP": "#FFC107",
        "MethylDackel": "#f0700e",
        "Modkit": "#B42CEA",
        "Pb-CpG-tools": "#8D9279",
        "Varlociraptor α = 0.1": "#004D40",
        "Varlociraptor α = 0.05": "#126e5f",
        "Varlociraptor α = 0.01": "#05AA8F",
        "Varlociraptor": "#05AA8F",
    }
    # Map nicer names
    df = df.copy()
    df["tool_name"] = df[color].map(meth_caller_to_name)
    charts = []
    df = df[df["platform"] != "Simulate"]
    for platform in df["platform"].unique():
        subset = df[df["platform"] == platform]
        # tools that appear in this subplot
        present_tools = subset["tool_name"].unique().tolist()
        color_domain = present_tools
        color_range = [tool_colors[t] for t in present_tools]
        ticks = list(np.logspace(0, np.log10(subset[x].max()), num=5))
        base = (
            alt.Chart(subset)
            .encode(
                x=alt.X(
                    f"{x}:Q",
                    title=x_title,
                    scale=alt.Scale(
                        type="symlog",
                        domain=[0, subset[x].max() * 1.1],
                    ),
                    axis=alt.Axis(values=ticks, labelAngle=-45),
                ),
                y=alt.Y(f"{y}:Q", title=y_title, scale=alt.Scale(type="log")),
                color=alt.Color(
                    "tool_name:N",
                    title="Caller",
                    scale=alt.Scale(domain=color_domain, range=color_range),
                ),
                shape=alt.Shape(f"{shape}:N", title="Task")
                if shape
                else alt.value("circle"),
                tooltip=(
                    [f"{x}:Q", f"{y}:Q", "tool_name:N", f"{shape}:N", "replicate:N"]
                ),
            )
            .mark_point(filled=False, size=30)
            .properties(height=height, width=height, title=platform)
            .interactive()
        )

        charts.append(base)

    long_reads = alt.hconcat(charts[1], charts[2]).resolve_scale(
        color="shared", y="shared"
    )
    return alt.hconcat(charts[0], long_reads).resolve_scale(color="independent")


# Read benchmark files from Snakemake input directory
records = []
benchmark_path = snakemake.input.benchmarks

# Validate that the benchmark path exists
if not os.path.exists(benchmark_path):
    raise FileNotFoundError(f"Benchmark directory not found: {benchmark_path}")

for root, _, files in os.walk(benchmark_path):
    for fname in files:
        if not fname.endswith(".bwa.benchmark.txt"):
            continue

        full_path = os.path.join(root, fname)
        try:
            df = pd.read_csv(full_path, sep="\t", usecols=["s", "max_rss"])
            p = Path(full_path)
            # Extract info from folder structure
            df["platform"] = p.parts[-4]  # → "Illumina_pe"
            df["meth_caller"] = p.parts[-3]  # → "varlociraptor"
            df["task"] = p.parts[-2]  # → "simulated_data_1"
            replicate = re.sub(
                r"_\d+-of-\d+", "", fname.replace(".bwa.benchmark.txt", "")
            )
            df["replicate"] = replicate
            # Clean up sample/replicate name
            # Compute avg of s and max_rss
            # Compute sum of s and max of max_rss
            df = df.groupby(
                ["platform", "meth_caller", "task", "replicate"], as_index=False
            ).agg(s=("s", "mean"), max_rss=("max_rss", "max"))

            records.append(df)
        except Exception as e:
            print(
                f"Warning: Failed to read benchmark file {full_path}: {e}",
                file=sys.stderr,
            )
if not records:
    raise ValueError(f"No benchmark files found in {benchmark_path}")

df_all = pd.concat(records, ignore_index=True)
df_all["platform"] = df_all["platform"].replace("Illumina_pe", "Illumina")

# df_all["task"] = np.where(
#     df_all["meth_caller"] != "varlociraptor",
#     "calling",
#     df_all["task"],
# )

# Define task to task_group mapping
task_group_mapping = {
    "bismark_align": "preprocessing",
    "deduplicate_bismark": "preprocessing",
    "bismark_methylation_extractor": "calling",
    "samtools_sort": "preprocessing",
    "samtools_merge": "preprocessing",
    "bissnp_compute": "calling",
    "bsmap_compute": "calling",
    "bsmap_extract": "calling",
    "methylDackel_compute_meth": "calling",
    "modkit": "calling",
    "pb-CpG-tools": "calling",
    "preprocessing": "preprocessing",
    "calling": "calling",
    "bissnp_extract": "calling",
}

# Assign task_group based on task name
df_all["task_group"] = df_all["task"].map(task_group_mapping).fillna("preprocessing")
# Compare different methylation calling callers
df_compare_callers = (
    df_all.groupby(
        ["platform", "meth_caller", "replicate", "task_group"], as_index=False
    )
    .agg(s=("s", "sum"), max_rss=("max_rss", "max"))
    .query("platform != 'multi_sample'")
    .assign(
        caller=lambda x: x["meth_caller"],
        minutes=lambda x: x["s"] / 60,
        max_rss_gb=lambda x: x["max_rss"] / 1024,
    )
    .groupby(["platform", "caller", "replicate", "task_group"], as_index=False)
    .agg({"minutes": "sum", "max_rss_gb": "max"})
    .assign(platform_caller=lambda x: x["platform"] + " - " + x["caller"])
)


# Create runtime and memory plots for all callers
runtime_chart = point_plot(
    df_compare_callers,
    x="minutes",
    y="max_rss_gb",
    color="caller",
    shape="task_group",
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
