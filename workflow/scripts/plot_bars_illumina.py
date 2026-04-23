import pickle
import sys
from hashlib import sha1

import altair as alt
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
alt.data_transformers.enable("vegafusion")


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
        distances_filtered = distances[(distances["sample"] == s) & (distances["meth_caller"] == m)]
        print(distances_filtered)
        mape_distance = (
            distances_filtered["mape"].values[0] if not distances_filtered.empty else 0.0
        )
        mae_distance = (
            distances_filtered["mae"].values[0] if not distances_filtered.empty else 0.0
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
                "distance": float(mae_distance),
                "distance_type": "Dₐ",
            }
        )
df_summary = pd.DataFrame(results)
colorblind_safe_palette = [
    "#D81B60",
    "#1E88E5",
    "#FFC107",
    "#05AA8F",
    "#004D40",
]
print(df_summary)
base = alt.Chart(df_summary).encode(
    x=alt.X("sample:N", axis=alt.Axis(labelAngle=-30), title=None),
    xOffset=alt.XOffset("meth_caller:N", sort=meth_callers),
    color=alt.Color(
        "meth_caller:N",
        title="Methylation caller",
        scale=alt.Scale(range=colorblind_safe_palette),
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
