import polars as pl
import altair as alt
import numpy as np
import pandas as pd

pd.set_option("display.max_rows", 1000)


truth_df = pl.read_csv(snakemake.input.truth)
replicate_df = pl.read_parquet(snakemake.input.results_rep)
meth_callers = snakemake.params.meth_callers


df = truth_df.join(
    replicate_df,
    left_on=["chrom", "pos"],
    right_on=["chromosome", "position"],
    how="inner",
).select(
    [
        pl.col("chrom"),
        pl.col("pos"),
        pl.col("methylation_level").alias("true_methylation"),
    ]
    + [pl.col(f"{caller}_methylation") for caller in meth_callers]
)


def compute_mape(true: pl.Series, pred: pl.Series) -> float:
    # Nulls entfernen (wichtig!)
    mask = pred.is_not_null()
    true = true.filter(mask)
    pred = pred.filter(mask)

    denom = np.maximum(true.to_numpy(), pred.to_numpy())
    denom[denom == 0] = 1.0

    mape = np.mean(np.abs(true.to_numpy() - pred.to_numpy()) / denom) * 100
    return float(mape)


mapes = {}
for caller in meth_callers:
    mape = compute_mape(df["true_methylation"], df[f"{caller}_methylation"])
    mapes[f"{caller}"] = mape
    print(f"{caller} MAPE: {mape}%")

long_df = df.melt(
    id_vars=["chrom", "pos", "true_methylation"],
    value_vars=[f"{caller}_methylation" for caller in meth_callers],
    variable_name="caller",
    value_name="caller_methylation",
)
bin_size = snakemake.params["bin_size"]

long_df = (
    long_df.with_columns(
        ((pl.col("true_methylation") // bin_size) * bin_size).alias("true_bin")
    )
    .with_columns(
        ((pl.col("caller_methylation") // bin_size) * bin_size).alias("caller_bin")
    )
    .group_by(["caller", "true_bin", "caller_bin"])
    .agg(pl.count().alias("count"))
)


# Polars → Pandas für Altair
heatmap_data = long_df.to_pandas()
# Alle Bins definieren
bins = np.arange(0, 101, bin_size)

# Alle Caller
callers = long_df["caller"].unique()
# Remove suffix "_methylation"
# Erstelle ein Grid aller Kombinationen
grid = pd.MultiIndex.from_product(
    [callers, bins, bins], names=["caller", "true_bin", "caller_bin"]
).to_frame(index=False)

# Merge mit den vorhandenen Counts, fehlende auf 0 setzen
heatmap_data_full = pd.merge(
    grid, heatmap_data, on=["caller", "true_bin", "caller_bin"], how="left"
).fillna({"count": 0})

max_count = heatmap_data_full["count"].max()

print(mapes)
print(heatmap_data_full.head())


def plot_heatmap(meth_caller, df):
    """Log-scaled heatmap for replicate methylation counts."""
    df = df[df["caller"] == f"{meth_caller}_methylation"]
    heatmap = (
        alt.Chart(
            df,
            title=alt.Title(
                f"{meth_caller}",
                subtitle=f" N = {df['count'].sum():.0f} MAPE = {mapes[meth_caller]:.2f}%",
            ),
        )
        .mark_rect()
        .encode(
            x=alt.X(
                "true_bin:O", sort=list(range(0, 101, bin_size)), title="Truth bins"
            ),
            y=alt.Y(
                "caller_bin:O",
                sort=list(range(100, -1, -bin_size)),
                title=f"{meth_caller} bins",
            ),
            color=alt.Color(
                "count:Q",
                scale=alt.Scale(type="log", scheme="viridis", domain=[1, max_count]),
            ),
            tooltip=["true_bin:O", "caller_bin:O", "count:Q"],
        )
    )
    return heatmap


heatmaps = []
callers = [c.replace("_methylation", "") for c in callers]

for meth_caller in callers:
    heatmap = plot_heatmap(meth_caller, heatmap_data_full)
    heatmaps.append(heatmap)
heatmap = alt.hconcat(*heatmaps)
# =============================
# 5️⃣ Speichern
# =============================

heatmap.save(snakemake.output[0])
