import altair as alt
import polars as pl

pl.Config.set_tbl_rows(1000)
pl.Config.set_tbl_cols(100)


def parse_cov(file_path, df):
    """Parses the coverage file and extracts coverage information."""
    cov_df = (
        pl.read_csv(
            file_path,
            separator="\t",
            has_header=False,
            schema={
                "chrom": pl.Utf8,
                "pos": pl.Int64,
                "end": pl.Int64,
                "coverage": pl.Float64,
            },
        ).select(["chrom", "pos", "coverage"])
    ).with_columns(
        pl.col("pos") + 1  # Convert to 1-based position
    )
    return df.join(
        cov_df,
        on=["chrom", "pos"],
        how="left",
    )


def plot_cov_dist(df_caller, caller, plot_type):
    chart = (
        alt.Chart(df_caller)
        .mark_bar()
        .encode(
            x=alt.X("coverage", bin=alt.Bin(step=5)),
            y=alt.Y("count()", stack=False),
        )
        .properties(title=f"{caller}: {plot_type}")
    )
    return chart



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
    "no_bias": "No bias",
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
    "No bias": "#888888",
}


def compute_precision_recall(df, meth_callers, bin_size=5):

    df = df.with_columns(
        pl.when(pl.col("true_methylation") > 0)
        .then(1)
        .otherwise(0)
        .alias("truth_binary")
    )

    results = []
    cov_charts = []
    print(df.filter((pl.col("no_bias_methylation") == 0) & (pl.col("methylDackel_methylation") > 0)))
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

        chart_tn_cv = plot_cov_dist(
            df_caller.filter(
                (pl.col(f"{caller}_binary") == 0) & (pl.col("truth_binary") == 0)
            ),
            caller,
            "True Negatives",
        )
        chart_tp_cv = plot_cov_dist(
            df_caller.filter(
                (pl.col(f"{caller}_binary") == 1) & (pl.col("truth_binary") == 1)
            ),
            caller,
            "True Positives",
        )
        chart_fn_cv = plot_cov_dist(
            df_caller.filter(
                (pl.col(f"{caller}_binary") == 0) & (pl.col("truth_binary") == 1)
            ),
            caller,
            "False Negatives",
        )
        chart_fp_cv = plot_cov_dist(
            df_caller.filter(
                (pl.col(f"{caller}_binary") == 1) & (pl.col("truth_binary") == 0)
            ),
            caller,
            "False Positives",
        )
        cov_charts.append(
            alt.concat(chart_tn_cv, chart_tp_cv, chart_fn_cv, chart_fp_cv)
        )
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

    return pl.DataFrame(results), alt.vconcat(*cov_charts)


alt.data_transformers.enable("vegafusion")
meth_callers = snakemake.params.meth_callers + ["no_bias"]
truth_df = pl.read_csv(snakemake.input[0], schema_overrides={"chrom": pl.Utf8})
tools_df = pl.read_parquet(snakemake.input[1]).with_columns(
    pl.col("chromosome").cast(pl.Utf8)
)
no_bias_df = pl.read_parquet(snakemake.input["no_bias"]).with_columns(
    pl.col("chromosome").cast(pl.Utf8)
)
print("Tools", tools_df.filter(pl.col("position").is_in([5030346])).head(20))
print("No bias", no_bias_df.filter(pl.col("position").is_in([5030346])).head(20))

df = tools_df.join(
    no_bias_df,
    on=["chromosome", "position"],
    how="inner",
).rename(
    {"tool_methylation": "no_bias_methylation"}
)
print("DF", df.head(20))

df = truth_df.join(
    df,
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
print("Combined", df.head(20))

df = parse_cov(snakemake.input[2], df)


metrics_df, cov_charts = compute_precision_recall(df, meth_callers)
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
cov_charts.save(snakemake.output[1])
