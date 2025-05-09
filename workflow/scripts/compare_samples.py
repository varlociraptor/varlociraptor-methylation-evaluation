import pandas as pd
import altair as alt
import sys

# Redirect error output to log file
sys.stderr = open(snakemake.log[0], "w")

# Show all columns during debugging
pd.set_option("display.max_columns", None)


def plot_correlation(df):
    num_samples = df["Sample"].nunique()
    num_comparisons = df["Comparison"].nunique()
    chart_width = num_samples * 100
    chart_height = num_comparisons * 30

    # Create heatmap
    heatmap = (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X("Sample:N", title="Sample"),
            y=alt.Y("Comparison:N", title="Comparison"),
            color=alt.Color(
                "Correlation:Q",
                scale=alt.Scale(scheme="redyellowgreen", domain=[0, 1]),
                legend=alt.Legend(title="Spearman ρ"),
            ),
            tooltip=[
                "Sample",
                "Comparison",
                alt.Tooltip("Correlation:Q", format=".4f"),
            ],
        )
        .properties(
            title="Spearman Correlation",
            width=chart_width,
            height=chart_height,
        )
    )

    # Overlay correlation values as text
    text = (
        alt.Chart(df)
        .mark_text(baseline="middle", fontSize=12)
        .encode(
            x="Sample:N",
            y="Comparison:N",
            text=alt.Text("Correlation:Q", format=".4f"),
            color=alt.value("white"),
        )
    )

    # Save chart
    return heatmap + text


def spearman_correlation_method(df, sample_name):
    """
    Calculates Spearman correlations between 'varlo_methylation' and
    all other *_methylation columns for a single sample.
    """

    correlations = []
    for method in snakemake.params["methods"]:
        if method == "varlo":
            continue
        corr = df["varlo_methylation"].corr(
            df[f"{method}_methylation"], method="spearman"
        )
        ref_tool = f"varlo vs. {method}"
        correlations.append(
            {
                "Sample": sample_name,
                "Comparison": ref_tool,
                "Correlation": corr,
            }
        )

    return pd.DataFrame(correlations)


def spearman_correlation_sample(df, sample_name):
    """
    Calculates Spearman correlations between 'varlo_methylation' and
    all other *_methylation columns for a single sample.
    """
    correlations = []
    for method in snakemake.params["methods"]:
        df_short = df.dropna(
            subset=[f"{method}_methylation_rep1", f"{method}_methylation_rep2"]
        )
        corr = df_short[f"{method}_methylation_rep1"].corr(
            df_short[f"{method}_methylation_rep2"], method="spearman"
        )
        correlations.append(
            {
                "Sample": sample_name,
                "Comparison": method,
                "Correlation": corr,
            }
        )

    return pd.DataFrame(correlations)


def plot_scatter_replicates(df_dict):
    scatter_plots = []

    for sample_name, df in df_dict.items():
        for method in snakemake.params["methods"]:
            x_col = f"{method}_methylation_rep1"
            y_col = f"{method}_methylation_rep2"

            if x_col in df.columns and y_col in df.columns:
                scatter = (
                    alt.Chart(df)
                    .mark_circle(size=10, opacity=0.4)
                    .encode(
                        x=alt.X(x_col, title=f"{method} Rep1"),
                        y=alt.Y(y_col, title=f"{method} Rep2"),
                        tooltip=["chromosome", "position", x_col, y_col],
                    )
                    .properties(
                        title=f"{sample_name} - {method}",
                        width=300,
                        height=300,
                    )
                )
                scatter_plots.append(scatter)

    return alt.hconcat(*scatter_plots)


# Read all sample files and compute correlations
replicate_dfs = {}
method_correlation_dfs = []
sample_correlation_dfs = []

for sample_file in snakemake.input["samples"]:
    replicate_name = sample_file.split("/")[-3]  # Extract protocol name as sample name
    df = pd.read_parquet(sample_file, engine="pyarrow")
    method_correlation_dfs.append(spearman_correlation_method(df, replicate_name))

    # Combine replicate DataFrames
    samplename = "_".join(replicate_name.split("_")[:-1])
    if samplename in replicate_dfs:
        replicate_dfs[samplename] = pd.merge(
            replicate_dfs[samplename],
            df,
            on=["chromosome", "position"],
            how="outer",
            suffixes=("_rep1", "_rep2"),
        )
    else:
        replicate_dfs[samplename] = df


for sample_name, df in replicate_dfs.items():
    # Calculate Spearman correlations for each sample
    # method_correlation_dfs.append(spearman_correlation_method(df, sample_name))
    sample_correlation_dfs.append(spearman_correlation_sample(df, sample_name))


# Combine all correlation DataFrames
correlation_all_methods = pd.concat(method_correlation_dfs, ignore_index=True)
correlation_all_samples = pd.concat(sample_correlation_dfs, ignore_index=True)
# print(correlation_all)

chart_methods = plot_correlation(correlation_all_methods)
chart_samples = plot_correlation(correlation_all_samples)
scatter_all = plot_scatter_replicates(replicate_dfs)
print(correlation_all_samples)
# Save charts
alt.vconcat(chart_methods, chart_samples, scatter_all).save(
    snakemake.output["plot"], scale_factor=2.0
)
