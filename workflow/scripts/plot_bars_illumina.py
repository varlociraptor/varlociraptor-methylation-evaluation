import pandas as pd
import altair as alt
import sys
import numpy as np
import pickle

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
alt.data_transformers.enable("vegafusion")


# -----------------------------
# Main execution
# -----------------------------
df = pd.read_parquet(snakemake.input["df"], engine="pyarrow")
mapes = pd.read_parquet(snakemake.input["mapes"], engine="pyarrow")
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

        # Filter mapes für Sample und Meth Caller
        mapes_filtered = mapes[(mapes["sample"] == s) & (mapes["meth_caller"] == m)]
        distance = mapes_filtered["mape"].values[0] if not mapes_filtered.empty else 0.0

        results.append(
            {
                "sample": s,
                "meth_caller": meth_caller_to_name.get(m, m),
                "number": str(number)[:3],
                "distance": float(distance),
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
bars = (
    alt.Chart(df_summary)
    .mark_bar()
    .encode(
        x=alt.X(
            "sample:N",
            axis=alt.Axis(labelAngle=-30),
            title=None,
        ),
        xOffset=alt.XOffset("meth_caller:N", sort=meth_callers),
        y=alt.Y("distance:Q", title="Discordance"),
        color=alt.Color(
            "meth_caller:N",
            title="Methylation caller",
            scale=alt.Scale(range=colorblind_safe_palette),
            sort=meth_callers,
        ),
        tooltip=["sample:N", "meth_caller:N", "distance:Q", "number:Q"],
    )
    .interactive()
)

labels = (
    alt.Chart(df_summary)
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

illumina_histo = bars + labels
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
