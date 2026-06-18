import altair as alt
import polars as pl

meth_caller_to_name = {
    "varlo_0.1": "Varlociraptor α = 0.1",
    "varlo_1.0": "Varlociraptor α = 1.0",
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
    "Varlociraptor α = 1.0": "#012b24",
    "Varlociraptor α = 0.1": "#004D40",
    "Varlociraptor α = 0.05": "#126e5f",
    "Varlociraptor α = 0.01": "#05AA8F",
}


def compute_precision_recall(df, meth_callers):
    df = df.with_columns(
        pl.when(pl.col("true_methylation") > 0)
        .then(1)
        .otherwise(0)
        .alias("truth_binary")
    )

    results = []

    for caller in meth_callers:
        df = df.with_columns(
            pl.when(pl.col(f"{caller}_methylation") > 0)
            .then(1)
            .otherwise(0)
            .alias(f"{caller}_binary")
        )

        TP = ((df[f"{caller}_binary"] == 1) & (df["truth_binary"] == 1)).sum()
        FP = ((df[f"{caller}_binary"] == 1) & (df["truth_binary"] == 0)).sum()
        FN = ((df[f"{caller}_binary"] == 0) & (df["truth_binary"] == 1)).sum()

        precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
        recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0

        results.append(
            {
                "caller": caller,
                "precision": precision,
                "recall": recall,
            }
        )

    return pl.DataFrame(results)


alt.data_transformers.enable("vegafusion")
meth_callers = snakemake.params.meth_callers
truth_df = pl.read_csv(snakemake.input[0])
tools_df = pl.read_parquet(snakemake.input[1])

df = truth_df.join(
    tools_df,
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


metrics_df = compute_precision_recall(df, meth_callers)
metrics_df = metrics_df.with_columns(
    pl.col("caller").replace(meth_caller_to_name).alias("tool_name")
)
print(metrics_df)
color_domain = (
    metrics_df.select(pl.col("tool_name").unique()).to_series().sort().to_list()
)
print(color_domain)
color_range = [tool_colors[t] for t in color_domain]
print(color_range)
chart = (
    alt.Chart(metrics_df.to_pandas())
    .mark_point(size=150)
    .encode(
        x=alt.X("recall:Q", title="Recall"),
        y=alt.Y("precision:Q", title="Precision"),
        color=alt.Color(
            "tool_name:N",
            scale=alt.Scale(
                domain=color_domain,
                range=color_range,
            ),
            title="Methylation caller",
        ),
        tooltip=["caller", "precision", "recall"],
    )
    .properties(width=500, height=400, title="Precision vs Recall")
)

chart.save(snakemake.output[0])
