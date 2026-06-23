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
        df_caller = df.filter(pl.col(f"{caller}_methylation").is_not_null())
        df_caller = df_caller.with_columns(
            pl.when(pl.col(f"{caller}_methylation") > 0)
            .then(1)
            .otherwise(0)
            .alias(f"{caller}_binary")
        )

        TP = (
            (df_caller[f"{caller}_binary"] == 1) & (df_caller["truth_binary"] == 1)
        ).sum()
        FP = (
            (df_caller[f"{caller}_binary"] == 1) & (df_caller["truth_binary"] == 0)
        ).sum()
        FN = (
            (df_caller[f"{caller}_binary"] == 0) & (df_caller["truth_binary"] == 1)
        ).sum()
        TN = (
            (df_caller[f"{caller}_binary"] == 0) & (df_caller["truth_binary"] == 0)
        ).sum()
        P = df_caller[f"{caller}_binary"].sum()
        N = df_caller.height - P
        precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
        recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
        positive_rate = TP / P if P > 0 else 0.0
        negative_rate = TN / N if N > 0 else 0.0

        print(
            f"caller: {caller}:",
            "\n\tTP: ",
            TP,
            "\n\tTN: ",
            TN,
            "\n\tP: ",
            P,
            "\n\tN: ",
            N,
            "\n\tpositive_rate: ",
            positive_rate,
            "\n\tnegative_rate: ",
            negative_rate,
        )
        results.append(
            {
                "caller": caller,
                "positive_rate": positive_rate,
                "negative_rate": negative_rate,
            }
        )

    return pl.DataFrame(results)


alt.data_transformers.enable("vegafusion")
meth_callers = snakemake.params.meth_callers
truth_df = pl.read_csv(snakemake.input[0], schema_overrides={"chrom": pl.Utf8})
tools_df = pl.read_parquet(snakemake.input[1]).with_columns(
    pl.col("chromosome").cast(pl.Utf8)
)
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
color_domain = (
    metrics_df.select(pl.col("tool_name").unique()).to_series().sort().to_list()
)
color_range = [tool_colors[t] for t in color_domain]
chart = (
    alt.Chart(metrics_df.to_pandas())
    .mark_point(size=150)
    .encode(
        x=alt.X("positive_rate:Q", title="Positive Rate"),
        y=alt.Y("negative_rate:Q", title="Negative Rate"),
        color=alt.Color(
            "tool_name:N",
            scale=alt.Scale(
                domain=color_domain,
                range=color_range,
            ),
            title="Methylation caller",
        ),
        tooltip=["caller", "positive_rate", "negative_rate"],
    )
    .properties(width=500, height=400, title="Positive Rate vs Negative Rate")
)

chart.save(snakemake.output[0])
