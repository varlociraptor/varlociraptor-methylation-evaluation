import pandas as pd
import altair as alt
import sys
import numpy as np

# Redirect error output to log file
sys.stderr = open(snakemake.log[0], "w")

# Show all columns during debugging
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)


def correlation_method(df, sample_name):
    """
    Calculates correlations between 'varlo_methylation' and
    all other *_methylation columns for a single sample.
    """

    correlations = []
    for method in snakemake.params["methods"]:
        if method == "varlo":
            continue
        pearson = df["varlo_methylation"].corr(
            df[f"{method}_methylation"], method="pearson"
        )
        spearman = df["varlo_methylation"].corr(
            df[f"{method}_methylation"], method="spearman"
        )
        ref_tool = f"varlo vs. {method}"
        correlations.append(
            {
                "sample": sample_name,
                "comparison": ref_tool,
                "pearson_corr": pearson,
                "spearman_corr": spearman,
            }
        )
    return pd.DataFrame(correlations)


def correlation_sample(df, sample_name):
    """
    Calculates correlations between 'varlo_methylation' and
    all other *_methylation columns for a single sample.
    """
    correlations = []
    for method in snakemake.params["methods"]:
        df_short = df.dropna(
            subset=[f"{method}_methylation_rep1", f"{method}_methylation_rep2"]
        )
        pearson = df_short[f"{method}_methylation_rep1"].corr(
            df_short[f"{method}_methylation_rep2"], method="pearson"
        )
        spearman = df_short[f"{method}_methylation_rep1"].corr(
            df_short[f"{method}_methylation_rep2"], method="spearman"
        )

        squared_errors = (
            df_short[f"{method}_methylation_rep1"]
            - df_short[f"{method}_methylation_rep2"]
        ) ** 2
        mean_squared_error = squared_errors.mean()
        rmse = np.sqrt(mean_squared_error)

        correlations.append(
            {
                "sample": sample_name,
                "comparison": method,
                "pearson_corr": pearson,
                "spearman_corr": spearman,
                "rmse_corr": rmse,
            }
        )

    return pd.DataFrame(correlations)


def plot_correlation(df, corr_method):
    num_samples = df["sample"].nunique()
    num_comparisons = df["comparison"].nunique()
    chart_width = num_samples * 100
    chart_height = num_comparisons * 30
    if corr_method == "rmse":
        domain = [10, 0]
    else:
        domain = [0, 1]
    # Create heatmap
    heatmap = (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X("sample:N", title="sample"),
            y=alt.Y("comparison:N", title="comparison"),
            color=alt.Color(
                f"{corr_method}_corr:Q",
                scale=alt.Scale(scheme="redyellowgreen", domain=[0, 1]),
                legend=alt.Legend(title=corr_method),
            ),
            # tooltip=[
            #     "sample",
            #     "comparison",
            #     alt.Tooltip(f"{corr_method}_corr:Q", format=".4f"),
            # ],
        )
        .properties(
            title=f"{corr_method.capitalize()} Correlation",
            width=chart_width,
            height=chart_height,
        )
    )

    # Overlay correlation values as text
    text = (
        alt.Chart(df)
        .mark_text(baseline="middle", fontSize=12)
        .encode(
            x="sample:N",
            y="comparison:N",
            text=alt.Text(f"{corr_method}_corr:Q", format=".4f"),
            color=alt.value("white"),
        )
    )

    # Save chart
    return heatmap + text


def plot_scatter_replicates(df_dict, corr_method):
    charts_per_sample = []

    for sample_name, df in df_dict.items():
        sample_charts = []
        for method in snakemake.params["methods"]:
            x_col = f"{method}_methylation_rep1"
            y_col = f"{method}_methylation_rep2"

            # if corr_method == "spearman":
            #     df[f"{method}_rank_rep1"] = df[x_col].rank()
            #     df[f"{method}_rank_rep2"] = df[y_col].rank()
            #     x_col = f"{method}_rank_rep1"
            #     y_col = f"{method}_rank_rep2"
            df_temp = df.dropna(subset=[x_col, y_col])
            # if x_col in df.columns and y_col in df.columns:
            scatter = (
                alt.Chart(df_temp)
                .mark_circle(size=10, opacity=0.4)
                .encode(
                    x=alt.X(x_col, title=f"{method} Rep1"),
                    y=alt.Y(y_col, title=f"{method} Rep2"),
                    # tooltip=["chromosome", "position", x_col, y_col],
                )
                .properties(
                    title=f"Datapoints {len(df_temp)}",
                    width=150,
                    height=150,
                )
            )
            sample_charts.append(scatter)

        if sample_charts:
            charts_per_sample.append(
                alt.hconcat(*sample_charts).properties(
                    title=f"{sample_name} {corr_method} Correlation",
                )
            )

    return alt.vconcat(*charts_per_sample)


# Read all sample files and compute correlations
replicate_dfs = {}
method_correlation_dfs = []
sample_correlation = []

for sample_file in snakemake.input["samples"]:
    replicate_name = sample_file.split("/")[-3]  # Extract protocol name as sample name
    df = pd.read_parquet(sample_file, engine="pyarrow")
    method_correlation_dfs.append(correlation_method(df, replicate_name))

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

    # Calculate correlations for each sample
    sample_correlation.append(correlation_sample(df, sample_name))


# Combine all correlation DataFrames
correlation_all_methods = pd.concat(method_correlation_dfs, ignore_index=True)
correlation_samples = pd.concat(sample_correlation, ignore_index=True)


chart_methods_pearson = plot_correlation(correlation_all_methods, "pearson")
chart_methods_spearman = plot_correlation(correlation_all_methods, "spearman")

chart_samples_pearson = plot_correlation(correlation_samples, "pearson")
chart_samples_spearman = plot_correlation(correlation_samples, "spearman")
chart_samples_rmse = plot_correlation(correlation_samples, "rmse")


with pd.HDFStore(snakemake.output["table"]) as store:
    for key, df in replicate_dfs.items():
        store[key] = df

# Horizontally combine method correlation charts
charts_methods = alt.hconcat(chart_methods_pearson, chart_methods_spearman).properties(
    title="Correlation of Varlo with reference tools",
)

# Horizontally combine sample correlation charts
charts_samples = alt.hconcat(
    chart_samples_pearson, chart_samples_spearman, chart_samples_rmse
).properties(
    title="Correlation of replicates",
)

# Vertically stack all charts
full_chart = alt.vconcat(
    charts_methods,
    charts_samples,
    # scatter_all_spearman,
)

# Save the full chart
full_chart.save(
    snakemake.output["plot"],
    embed_options={"actions": False},
    inline=False,  # <- wichtig
)
