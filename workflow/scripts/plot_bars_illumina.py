import pickle
import sys
from hashlib import sha1

import altair as alt
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
alt.data_transformers.enable("vegafusion")

meth_caller_to_name = {
    "varlo_0.1": "Varlociraptor α = 0.1",
    "varlo_0.05": "Varlociraptor α = 0.05",
    "varlo_0.01": "Varlociraptor α = 0.01",
    "varlociraptor": "Varlociraptor",
    "bismark": "Bismark",
    "bsMap": "BSMAPz",
    "methylDackel": "MethylDackel",
    "modkit": "Modkit",
    "pb_CpG_tools": "Pb-CpG-tools",
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
}

# -----------------------------
# Main execution
# -----------------------------
df = pd.read_parquet(snakemake.input["df"], engine="pyarrow")
distances = pd.read_parquet(snakemake.input["distances"], engine="pyarrow")
bin_size = snakemake.params["bin_size"]
meth_callers = df["meth_caller"].unique().tolist()
samples = df["sample"].unique().tolist()
plot_type = snakemake.params.get("plot_type")
meth_caller_to_name = {
    "bismark": "Bismark",
    "bsMap": "BSMAPz",
    "bisSNP": "BisSNP",
    "methylDackel": "MethylDackel",
}

for m in meth_callers:
    if m.startswith("varlo_"):
        alpha = m.split("_")[1]
        meth_caller_to_name[m] = f"Varlociraptor α = {alpha}"
results = []
for s in samples:
    for m in meth_callers:
        # Filter df für Sample und Meth Caller
        df_filtered = df[(df["sample"] == s) & (df["meth_caller"] == m)]
        number = df_filtered["count"].sum() if not df_filtered.empty else 0

        # Filter distances für Sample und Meth Caller
        distances_filtered = distances[
            (distances["sample"] == s) & (distances["meth_caller"] == m)
        ]
        mape_distance = (
            distances_filtered["mape"].values[0]
            if not distances_filtered.empty
            else 0.0
        )
        mae_distance = (
            distances_filtered["mae"].values[0] if not distances_filtered.empty else 0.0
        )
        binary_concordance = (
            distances_filtered["binary_concordance"].values[0]
            if not distances_filtered.empty
            else 0.0
        )

        results.append(
            {
                "sample": s,
                "meth_caller": meth_caller_to_name.get(m, m),
                "number": str(number)[:3],
                "distance": float(mape_distance),
                "distance_type": "Dᵣ",
            }
        )
        results.append(
            {
                "sample": s,
                "meth_caller": meth_caller_to_name.get(m, m),
                "number": str(number)[:3],
                "distance": float(binary_concordance),
                "distance_type": "Bc",
            }
        )
        results.append(
            {
                "sample": s,
                "meth_caller": meth_caller_to_name.get(m, m),
                "number": str(number)[:3],
                "distance": float(mae_distance),
                "distance_type": "Dₐ",
            }
        )
df_summary = pd.DataFrame(results)
df_summary["tool_name"] = df_summary["meth_caller"].replace(meth_caller_to_name)
color_domain = sorted(df_summary["tool_name"].unique())
color_range = [tool_colors[t] for t in color_domain]
base = alt.Chart(df_summary).encode(
    x=alt.X("sample:N", axis=alt.Axis(labelAngle=-30), title=None),
    xOffset=alt.XOffset("meth_caller:N", sort=meth_callers),
    color=alt.Color(
        "meth_caller:N",
        title="Methylation caller",
        scale=alt.Scale(
            domain=color_domain,
            range=color_range,
        ),
        sort=meth_callers,
    ),
    tooltip=["sample:N", "meth_caller:N", "distance:Q", "number:Q"],
)

# Dr (hinterer Balken)
bars_dr = (
    base.transform_filter(alt.datum.distance_type == "Dᵣ")
    .mark_bar(opacity=0.7)
    .encode(y=alt.Y("distance:Q", title="Discordance"))
)

# Da (vorderer Balken, schraffiert)
bars_da = (
    base.transform_filter(alt.datum.distance_type == "Dₐ")
    .mark_bar(
        stroke="black",
        strokeWidth=1,
        strokeDash=[4, 2],  # "Schraffur"-Ersatz
    )
    .encode(y="distance:Q")
)

labels = (
    alt.Chart(df_summary)
    .transform_filter(alt.datum.distance_type == "Dᵣ")
    .mark_text(size=8, dy=-5, color="black")
    .transform_calculate(text_k="datum.number + 'k'")
    .encode(
        text="text_k:N",
        x="sample:N",
        xOffset=alt.XOffset("meth_caller:N", sort=meth_callers),
        y="distance:Q",
        tooltip=["sample:N", "meth_caller:N", "distance:Q", "number:Q"],
    )
    .interactive()
)

illumina_histo = (bars_dr + bars_da + labels).interactive()
if plot_type == "parquet":
    df_summary.to_parquet(snakemake.output[0])
elif plot_type == "pkl":
    with open(snakemake.output[0], "wb") as f:
        pickle.dump(illumina_histo, f)
else:
    illumina_histo.save(
        snakemake.output[0],
        embed_options={"actions": False},
        inline=False,
    )
