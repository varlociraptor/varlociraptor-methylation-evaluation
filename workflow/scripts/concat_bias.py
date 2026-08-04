from turtle import color

import altair as alt
import pandas as pd
import polars as pl


pl.Config.set_tbl_cols(100)

def bias_plots(df_long: pd.DataFrame, fdr: str):
    """Create bias, AF, and DP plots from long-format data."""
    # Bias category plot
    bias_chart = (
        alt.Chart(df_long.to_pandas())
        .mark_bar()
        .encode(
            x=alt.X(
                "category:N",
                axis=alt.Axis(labelAngle=-45),
                title=None,
                scale=alt.Scale(
                    domain=["Bias both reps", "Bias, AF = 0", "Bias, AF > 0"],
                ),
            ),
            y="count():Q",
            color=alt.Color(
                "bias_type_label:N",
                scale=alt.Scale(
                    domain=df_long["bias_type_label"].unique(),
                    range=["#D81B60", "#1E88E5"],
                ),
                title="Bias Type",
                # legend=None if platform_label != "Nanopore" else alt.Legend(),
            ),
            tooltip=["category", "count()", "bias_type_label"],
            column=alt.Column("platform_label:N", title=None),
        )
    )

    return bias_chart


def depth_plots(df: pd.DataFrame, fdr: str):
    df = (
        df.filter(pl.col("category") == "Bias, AF > 0")
        .with_columns(pl.max_horizontal("AF_rep1", "AF_rep2").round(2).alias("AF"))
        .filter((pl.col("DP_rep1") <= 500) & (pl.col("DP_rep2") <= 500))
        .with_columns(
            pl.when(pl.col("AF_rep1") == 0).then(pl.col("DP_rep1")).otherwise(pl.col("DP_rep2")).alias("DP_bias"),
            pl.when(pl.col("AF_rep1") > 0).then(pl.col("DP_rep1")).otherwise(pl.col("DP_rep2")).alias("DP_AF"),
        )
    )
    print(df)
    depth_chart = (
        alt.Chart(df.to_pandas())
        .mark_circle()
        .encode(
            x=alt.X("DP_bias:Q", title="Depth associated with bias"),
            y=alt.Y("DP_AF:Q", title="Depth associated with AF > 0"),
            tooltip=["DP_bias:Q", "DP_AF:Q"],
            color=alt.Color(
                "AF:Q",
                title="AF of non-bias sample",
                scale=alt.Scale(scheme="viridis"),  # oder "turbo", "plasma", "inferno"
            ),
            column=alt.Column("platform_label:N", title=None),
        )
        .properties(title="Depth scatter plot for 'Bias, AF > 0'")

    )
    return depth_chart


df_illumina = pl.read_parquet(snakemake.input["illumina"])
df_pacbio = pl.read_parquet(snakemake.input["pacbio"])
df_nanopore = pl.read_parquet(snakemake.input["nanopore"])


df_illumina = df_illumina.with_columns(pl.lit("Illumina").alias("platform_label"))
df_pacbio = df_pacbio.with_columns(pl.lit("PacBio").alias("platform_label"))
df_nanopore = df_nanopore.with_columns(pl.lit("Nanopore").alias("platform_label"))


df = pl.concat([df_illumina, df_pacbio, df_nanopore])
# Concat and create column platform_label for the corresponding platform
bias_chart = bias_plots(df, snakemake.params.fdr)
bias_chart.save(snakemake.output[0])

depth_chart = depth_plots(df, snakemake.params.fdr)
depth_chart.save(snakemake.output[1])
