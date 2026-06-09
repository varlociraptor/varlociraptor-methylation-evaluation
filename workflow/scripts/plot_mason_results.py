import altair as alt
import numpy as np
import pandas as pd
import polars as pl

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_rows", 1000)
pl.Config.set_tbl_cols(100)


truth_df = pl.read_csv(
    snakemake.input.truth,
    schema_overrides={"chrom": pl.Utf8},
)
replicate_df = pl.read_parquet(snakemake.input.results_rep)
meth_callers = snakemake.params.meth_callers
print(truth_df.head(), file=sys.stderr)
print(replicate_df.head(), file=sys.stderr)
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

(print(df.head(), file=sys.stderr),)


def compute_mape(df, meth_caller) -> float:
    print(df.head(), file=sys.stderr)
    df = df.with_columns(
        pl.max_horizontal(
            pl.col(f"{meth_caller}_methylation"),
            pl.col("true_methylation"),
        ).alias("denom")
    ).with_columns(
        pl.when(pl.col("denom") == 0)
        .then(0.0)
        .otherwise(
            (pl.col(f"{meth_caller}_methylation") - pl.col("true_methylation")).abs()
            / pl.col("denom")
        )
        .alias("mape_row")
    )
    mape = df.select(pl.col("mape_row").mean() * 100).item()
    return float(mape)


def compute_mae(df, meth_caller) -> float:
    df = df.with_columns(
        (pl.col(f"{meth_caller}_methylation") - pl.col("true_methylation"))
        .abs()
        .alias("mae_row")
    )

    mae = df.select(pl.col("mae_row").mean()).item()
    return float(mae)


distance_rows = []

for caller in meth_callers:
    distance_rows.append(
        {
            "meth_caller": caller,
            "mape": compute_mape(df, caller),
            "mae": compute_mae(df, caller),
        }
    )

distance_df = pl.DataFrame(distance_rows)

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


heatmap_data = long_df.to_pandas()
bins = np.arange(0, 101, bin_size)

callers = long_df["caller"].unique()
# Remove suffix "_methylation"
# Grid for all combinations
grid = pd.MultiIndex.from_product(
    [callers, bins, bins], names=["caller", "true_bin", "caller_bin"]
).to_frame(index=False)

heatmap_data_full = pd.merge(
    grid, heatmap_data, on=["caller", "true_bin", "caller_bin"], how="left"
).fillna({"count": 0})

max_count = heatmap_data_full["count"].max()


def plot_heatmap(meth_caller, df, distance_df):
    """Log-scaled heatmap for replicate methylation counts."""
    df = df[df["caller"] == f"{meth_caller}_methylation"]
    heatmap = (
        alt.Chart(
            df,
            title=alt.Title(
                f"{meth_caller}",
                subtitle=f" N = {df['count'].sum():.0f} Dᵣ = {distance_df.filter(pl.col('meth_caller') == meth_caller)['mape'].item():.2f}%, Dₐ = {distance_df.filter(pl.col('meth_caller') == meth_caller)['mae'].item():.2f}%",
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
for meth_caller in meth_callers:
    heatmap = plot_heatmap(meth_caller, heatmap_data_full, distance_df)
    heatmaps.append(heatmap)
heatmap = alt.hconcat(*heatmaps)
print(heatmap_data_full)

heatmap_data_full["distance"] = (
    heatmap_data_full["caller_bin"] - heatmap_data_full["true_bin"]
).abs()

distance_plot_df = heatmap_data_full.groupby(["caller", "distance"], as_index=False)[
    "count"
].sum()
print(distance_plot_df)
distance_plot = (
    alt.Chart(distance_plot_df)
    .mark_line()
    .encode(
        x=alt.X(
            "distance:O",
            title="caller_bin − true_bin",
            sort="ascending",
        ),
        y=alt.Y(
            "count:Q",
            title="Count",
        ),
        tooltip=["distance:O", "count:Q"],
        color=alt.Color(
            "caller:O",
            scale=alt.Scale(scheme="category10"),
            title="Methylation Caller",
        ),
    )
)
plot = alt.vconcat(heatmap, distance_plot).resolve_scale(color="independent")

plot.save(snakemake.output[0])
