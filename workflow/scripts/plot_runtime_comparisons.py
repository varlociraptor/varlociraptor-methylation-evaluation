import os
import pandas as pd
import re
import altair as alt


all_data = []

path = snakemake.input.benchmarks

for root, dirs, files in os.walk(path):
    for f in files:
        if not f.endswith(".benchmark.txt"):
            continue

        full_path = os.path.join(root, f)
        df = pd.read_csv(full_path, sep="\t", usecols=["s", "max_rss"])

        # Plattform und Methode
        df["meth_caller"] = os.path.basename(os.path.dirname(full_path))
        df["platform"] = os.path.basename(os.path.dirname(os.path.dirname(full_path)))

        # Sample-Name: _x-of-y entfernen
        base_name = os.path.basename(f)
        base_name = re.sub(r"_\d+-of-\d+", "", base_name)
        sample_name = base_name.replace(".benchmark.txt", "")
        df["replicate"] = sample_name

        all_data.append(df)

# Alle Files zusammenführen
df_all = pd.concat(all_data, ignore_index=True)

# Zusammenfassen nach Platform, Methode und Sample
df_summary = df_all.groupby(
    ["platform", "meth_caller", "replicate"], as_index=False
).agg({"s": "sum", "max_rss": "max"})

print(df_summary.to_string())

# Tool zusammenfassen: calling + preprocessing → Varlociraptor
df_summary["tool"] = df_summary["meth_caller"].replace(
    {"calling": "Varlociraptor", "preprocessing": "Varlociraptor"}
)

# Aggregieren: s summieren, max_rss maximal
df_agg = df_summary.groupby(["platform", "tool", "replicate"], as_index=False).agg(
    {"s": "sum", "max_rss": "max"}
)

# Minuten und GB berechnen
df_agg["minutes"] = df_agg["s"] / 60
df_agg["max_rss_gb"] = df_agg["max_rss"] / 1024

# Neue x-Achse: Plattform + Tool
df_agg["platform_tool"] = df_agg["platform"] + " - " + df_agg["tool"]

# Tooltip: Substeps kann man weglassen, da aggregiert
df_agg["replicate_tool"] = df_agg["replicate"] + " (" + df_agg["tool"] + ")"

# Runtime Boxplot
runtime_plot = (
    alt.Chart(df_agg)
    .mark_boxplot(extent="min-max")
    .encode(
        x=alt.X("platform_tool:N", title="Platform - Tool"),
        y=alt.Y("minutes:Q", title="Runtime (min)", scale=alt.Scale(type="log")),
        color="platform:N",
    )
    .properties(width=600, height=400, title="Runtime by Platform and Tool")
)

runtime_points = (
    alt.Chart(df_agg)
    .mark_point(filled=True, size=50)
    .encode(
        x="platform_tool:N",
        y="minutes:Q",
        color="platform:N",
        tooltip=["replicate_tool", "minutes"],
    )
)

runtime_chart = runtime_plot + runtime_points

# Memory Boxplot
memory_plot = (
    alt.Chart(df_agg)
    .mark_boxplot(extent="min-max")
    .encode(
        x=alt.X("platform_tool:N", title="Platform - Tool"),
        y=alt.Y("max_rss_gb:Q", title="Max RSS (GB)", scale=alt.Scale(type="log")),
        color="platform:N",
    )
    .properties(width=600, height=400, title="Memory Usage by Platform and Tool")
)

memory_points = (
    alt.Chart(df_agg)
    .mark_point(filled=True, size=50)
    .encode(
        x="platform_tool:N",
        y="max_rss_gb:Q",
        color="platform:N",
        tooltip=["replicate_tool", "max_rss_gb"],
    )
)

memory_chart = memory_plot + memory_points

# Charts zusammenführen
chart = alt.hconcat(runtime_chart, memory_chart).resolve_scale(color="independent")

# Chart speichern
chart.save(
    snakemake.output[0],
    embed_options={"actions": False},
    inline=False,
)
