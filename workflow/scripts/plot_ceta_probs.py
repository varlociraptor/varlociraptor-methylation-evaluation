import pandas as pd
from pathlib import Path
from functools import reduce
import altair as alt
import polars as pl
import sys

# sys.stderr = open(snakemake.log[0], "w")

pl.Config.set_tbl_rows(-1)  # show all rows
pl.Config.set_tbl_cols(-1)  # show all columns

files = snakemake.input


method_to_name = {
    "MethylSeq_HG002_LAB01_REP01": "MethylSeq",
    "EMSeq_HG002_LAB01_REP01": "EMSeq",
    "ceta_multi": "EMSeq + MethylSeq",
    "untreated": "Untreated",
}

colorblind_safe_palette = [
    "#D81B60",
    "#1E88E5",
    "#FFC107",
    "#05AA8F",
    "#004D40",
]

dfs = []

for f in files:
    method = method_to_name.get(Path(f).parts[-3], "Unknown")

    # Load parquet
    df = pl.scan_parquet(f)
    df = df.with_columns(pl.lit(method).alias("method"))
    dfs.append(df)

df = pl.concat(dfs)
bin_size = 20
intervals = [i / bin_size for i in range(bin_size)]
interval_labels = ["missing"] + [
    f"{i / bin_size} - { (i + 1) / bin_size}" for i in range(bin_size)
]

df = df.with_columns(
    pl.col("prob_present")
    .cut(intervals, labels=interval_labels, left_closed=True)
    .alias("prob_present_bin"),
    pl.col("prob_absent")
    .cut(intervals, labels=interval_labels, left_closed=True)
    .alias("prob_absent_bin"),
    pl.col("prob_artifact")
    .cut(intervals, labels=interval_labels, left_closed=True)
    .alias("prob_artifact_bin"),
)



# print(
#     df.filter(pl.col("prob_absent_bin").is_in([32225333, 14102609]))
#       .collect()
#       .head(200)
# )

# df = df.with_columns(pl.all().fill_null("0"))

# ---- Funktion für die Histogramme ----
def make_plot(df, category):
    df = df.sort("method", f"prob_{category}_bin")
    # df = df.filter(pl.col(f"prob_{category}_bin").is_null())
    # print(
    print(df.collect().tail(20))
 
    
    df = df.collect().to_pandas()
    col = f"prob_{category}_bin"
    df[col] = df[col].cat.add_categories(["missing"]).fillna("missing")

    df[f"prob_{category}_bin"] = df[f"prob_{category}_bin"].fillna("missing")
    print(df.tail(20))
    # df.fillna(-1, inplace=True)
    
    # df = df.sort_values(by=["method"])
    chart = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X(
                f"prob_{category}_bin:N",
                title=None,
                scale=alt.Scale(domain=interval_labels),
                axis=alt.Axis(labelAngle=-45),
            ),
            xOffset=alt.XOffset("method:N"),
            y=alt.Y("count()", title="Count"),
            color=alt.Color("method:N", title="Method", scale=alt.Scale(range=colorblind_safe_palette)),
            tooltip=["method:N", f"prob_{category}_bin:N", "count()"],
        )
        .properties(title=f"Distribution of prob_{category}")
    )
    return chart


# ---- Drei Plots erzeugen ----
plot_present = make_plot(df, "present")
plot_absent = make_plot(df, "absent")
plot_artifact = make_plot(df, "artifact")

heatmap_plots = alt.vconcat(plot_present, plot_absent, plot_artifact).resolve_scale(
    y="shared"
)


heatmap_plots.save(snakemake.output[0], embed_options={"actions": False}, inline=False)
