import pandas as pd
import altair as alt
import sys

# Redirect error output to log file
sys.stderr = open(snakemake.log[0], "w")

# Show all columns during debugging
pd.set_option("display.max_columns", None)


def pearson_correlation_single(df, sample_name):
    """
    Calculates Pearson correlations between 'varlo_methylation' and
    all other *_methylation columns for a single sample.
    """
    methylation_cols = [
        col
        for col in df.columns
        if col.endswith("_methylation") and col != "varlo_methylation"
    ]

    correlations = []
    for col in methylation_cols:
        corr = df["varlo_methylation"].corr(df[col])
        ref_tool = f"varlo vs. {col.replace('_methylation', '')}"
        correlations.append(
            {
                "Sample": sample_name,
                "Comparison": ref_tool,
                "Correlation": corr,
            }
        )

    return pd.DataFrame(correlations)


# Read all sample files and compute correlations
correlation_dfs = []
for sample_file in snakemake.input["samples"]:
    sample_name = sample_file.split("/")[2]  # Extract protocol name as sample name
    df = pd.read_parquet(sample_file, engine="pyarrow")
    correlation_dfs.append(pearson_correlation_single(df, sample_name))

# Combine all correlation DataFrames
correlation_all = pd.concat(correlation_dfs, ignore_index=True)

num_samples = correlation_all["Sample"].nunique()
num_comparisons = correlation_all["Comparison"].nunique()
chart_width = num_samples * 100
chart_height = num_comparisons * 30

# Create heatmap
heatmap = (
    alt.Chart(correlation_all)
    .mark_rect()
    .encode(
        x=alt.X("Sample:N", title="Sample"),
        y=alt.Y("Comparison:N", title="Comparison"),
        color=alt.Color(
            "Correlation:Q",
            scale=alt.Scale(scheme="redyellowgreen", domain=[0, 1]),
            legend=alt.Legend(title="Pearson r"),
        ),
        tooltip=[
            "Sample",
            "Comparison",
            alt.Tooltip("Correlation:Q", format=".4f"),
        ],
    )
    .properties(
        title="Pearson Correlation per Sample",
        width=chart_width,
        height=chart_height,
    )
)

# Overlay correlation values as text
text = (
    alt.Chart(correlation_all)
    .mark_text(baseline="middle", fontSize=12)
    .encode(
        x="Sample:N",
        y="Comparison:N",
        text=alt.Text("Correlation:Q", format=".4f"),
        color=alt.value("white"),
    )
)

# Save chart
(heatmap + text).save(snakemake.output["plot"], scale_factor=2.0)
