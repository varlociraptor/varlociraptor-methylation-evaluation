import os
import re
import pandas as pd
import altair as alt

def boxplot_with_points(df, x, y, color, x_title, y_title, title):
    """Reusable Altair boxplot + scatter overlay."""
    base = alt.Chart(df).encode(
        x=alt.X(f"{x}:N", title=x_title),
        y=alt.Y(f"{y}:Q", title=y_title, scale=alt.Scale(type="log")),
        color=f"{color}:N"
    )

    box = base.mark_boxplot(extent="min-max")
    points = base.mark_point(filled=True, size=50).encode(tooltip=[x, y])

    return (box + points).properties(width=600, height=400, title=title)


# -----------------------
# 🧠 Collect benchmark data
# -----------------------

records = []
path = snakemake.input.benchmarks

for root, _, files in os.walk(path):
    for f in files:
        if not f.endswith(".benchmark.txt"):
            continue

        full_path = os.path.join(root, f)
        df = pd.read_csv(full_path, sep="\t", usecols=["s", "max_rss"])

        # Extract metadata
        df["meth_caller"] = os.path.basename(os.path.dirname(full_path))
        df["platform"] = os.path.basename(os.path.dirname(os.path.dirname(full_path)))

        # Clean sample name
        replicate = re.sub(r"_\d+-of-\d+", "", f.replace(".benchmark.txt", ""))
        df["replicate"] = replicate

        records.append(df)

df_all = pd.concat(records, ignore_index=True)

# -----------------------
# 🔧 Aggregate results
# -----------------------

df_summary = (
    df_all.groupby(["platform", "meth_caller", "replicate"], as_index=False)
    .agg(s=("s", "sum"), max_rss=("max_rss", "max"))
)

# -----------------------
# 🧩 Compare tools
# -----------------------

# Instead of replace() on Series, use np.where (faster and simpler)
df_compare_tools = (
    df_summary.query("platform != 'multi_sample'")
    .assign(
        tool=lambda x: x["meth_caller"].replace(
            {"calling": "Varlociraptor", "preprocessing": "Varlociraptor"}
        ),
        minutes=lambda x: x["s"] / 60,
        max_rss_gb=lambda x: x["max_rss"] / 1024,
    )
    .groupby(["platform", "tool", "replicate"], as_index=False) .agg({"minutes": "sum", "max_rss_gb": "max"})
    .assign(platform_tool=lambda x: x["platform"] + " - " + x["tool"])
)

df_compare_varlo = (
    df_summary.query("meth_caller in ['calling', 'preprocessing']")
    .assign(
        minutes=lambda x: x["s"] / 60,
        max_rss_gb=lambda x: x["max_rss"] / 1024,
        platform_tool=lambda x: x["platform"] + " - " + x["meth_caller"],
    )
)

# -----------------------
# 📊 Plot charts
# -----------------------
print(df_compare_tools.to_string())

runtime_chart = boxplot_with_points(
    df_compare_tools,
    x="platform_tool", y="minutes", color="platform",
    x_title="Platform - Tool", y_title="Runtime (min)",
    title="Runtime of methylation calling tools"
)

memory_chart = boxplot_with_points(
    df_compare_tools,
    x="platform_tool", y="max_rss_gb", color="platform",
    x_title="Platform - Tool", y_title="Max RSS (GB)",
    title="Memory usage of methylation calling tools"
)

runtime_varlo_chart = boxplot_with_points(
    df_compare_varlo,
    x="platform_tool", y="minutes", color="platform",
    x_title="Platform - Step", y_title="Runtime (min)",
    title="Runtime of Varlociraptor steps"
)

# memory_varlo_chart = boxplot_with_points(
#     df_compare_varlo,
#     x="platform_tool", y="max_rss_gb", color="platform",
#     x_title="Platform - Step", y_title="Max RSS (GB)",
#     title="Memory usage of Varlociraptor steps"
# )

# Combine charts
tool_chart = alt.hconcat(runtime_chart, memory_chart)

# Save
tool_chart.save(
    snakemake.output["tools"],
    embed_options={"actions": False},
    inline=False,
)


runtime_varlo_chart.save(
    snakemake.output["varlo"],
    embed_options={"actions": False},
    inline=False,
)

