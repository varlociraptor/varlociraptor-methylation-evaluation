import pandas as pd
from pathlib import Path
from functools import reduce
import altair as alt


files = snakemake.input

method_to_name = {
    "MethylSeq_HG002_LAB01_REP01": "MethylSeq",
    "EMSeq_HG002_LAB01_REP01": "EMSeq",
    "ceta_multi": "EMSeq + MethylSeq",
    "untreated": "Untreated",
}
dfs = []

for f in files:
    method = method_to_name.get(Path(f).parts[-3], "Unknown")

    # Load parquet
    df = pd.read_parquet(f)
    df["method"] = method

    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

print(df.head())


# ---- Funktion für die Histogramme ----
def make_plot(category):

    # Compute normalized frequencies per method
    df_hist = (
        df.assign(bin=lambda d: pd.cut(d[f"prob_{category}"], bins=10))
        .groupby(["method", "bin"])
        .size()
        .groupby(level=0)  # normalize per method
        .apply(lambda x: x / x.sum())
    )

    print(df_hist.to_string())

    chart = (
        alt.Chart(df_hist)
        .mark_bar()
        .encode(
            x=alt.X("bin:N", title=f"{category} probability bin", sort=None),
            y=alt.Y("rel_freq:Q", title="Relative frequency"),
            color=alt.Color("method:N", title="Method"),
        )
        .properties(width=200, height=250, title=f"Distribution of prob_{category}")
    )
    return chart


# ---- Drei Plots erzeugen ----
plot_present = make_plot("present")
plot_absent = make_plot("absent")
plot_artifact = make_plot("artifact")

heatmap_plots = alt.hconcat(plot_present, plot_absent, plot_artifact).resolve_scale(
    x="independent", y="independent", color="independent"
)


heatmap_plots.save(snakemake.output[0], embed_options={"actions": False}, inline=False)
