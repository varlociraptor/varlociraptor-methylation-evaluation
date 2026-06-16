import altair as alt
import pandas as pd


def bias_plots(df_long: pd.DataFrame, fdr: str):
    """Create bias, AF, and DP plots from long-format data."""
    # Bias category plot
    bias_chart = (
        alt.Chart(df_long)
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


df_illumina = pd.read_parquet(snakemake.input["illumina"])
df_pacbio = pd.read_parquet(snakemake.input["pacbio"])
df_nanopore = pd.read_parquet(snakemake.input["nanopore"])

df_illumina["platform_label"] = "Illumina"
df_pacbio["platform_label"] = "PacBio"
df_nanopore["platform_label"] = "Nanopore"

df = pd.concat([df_illumina, df_pacbio, df_nanopore], ignore_index=True)
# Concat and create column platform_label for the corresponding platform
bias_chart = bias_plots(df, snakemake.params.fdr)
bias_chart.save(snakemake.output[0])
