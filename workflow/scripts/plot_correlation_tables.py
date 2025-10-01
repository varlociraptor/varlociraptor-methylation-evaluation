import pandas as pd
import altair as alt
import sys
import numpy as np
import re
from typing import Literal

# Redirect error output to log file
sys.stderr = open(snakemake.log[0], "w")

# Show all columns during debugging
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)


from typing import Literal


# https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
def compute_ccc(x, y):

    mean_x = x.mean()
    mean_y = y.mean()

    var_x = x.var(ddof=0)
    var_y = y.var(ddof=0)

    pearson = x.corr(y, method="pearson")

    ccc = (2 * pearson * (var_x**0.5) * (var_y**0.5)) / (
        var_x + var_y + (mean_x - mean_y) ** 2
    )

    return ccc


def compute_correlation(
    df: pd.DataFrame,
    sample_name: str,
    mode: Literal["meth_caller", "replicate"] = "meth_caller",
    corr_methods: list[str] = ["rmse"],
) -> pd.DataFrame:
    """
    Compute correlations and RMSE for methylation data.
    mode:
    - meth_caller: Compare the results of varlociraptor against the other tools on one replicate
    - replicate: Compare the results of a single tool on two replicates
    """
    correlations = []

    for meth_caller in snakemake.params["meth_callers"]:
        if mode == "meth_caller":
            if meth_caller == "varlo":
                continue
            col1 = "varlo_methylation"
            col2 = f"{meth_caller}_methylation"
            comp_name = f"varlo vs. {meth_caller}"

        elif mode == "replicate":
            col1 = f"{meth_caller}_methylation_rep1"
            col2 = f"{meth_caller}_methylation_rep2"
            if col1 not in df.columns or col2 not in df.columns:
                continue
            comp_name = meth_caller

        else:
            raise ValueError(f"Invalid mode: {mode}")

        # Drop NA rows for fair comparison
        df = df.dropna(subset=[col1, col2])
        if df.empty:
            continue

        metric_funcs = {
            "pearson": lambda: df[col1].corr(df[col2], method="pearson"),
            "spearman": lambda: df[col1].corr(df[col2], method="spearman"),
            "rmse": lambda: np.sqrt(((df[col1] - df[col2]) ** 2).mean()) / 100.0,
            "mape": lambda: (
                np.where(
                    (df[col1] == 0) & (df[col2] == 0),
                    0,
                    np.abs(df[col1] - df[col2]) / df[[col1, col2]].max(axis=1),
                ).mean()
                * 100.0
            ),
            "ccc": lambda: compute_ccc(df[col1], df[col2]),
        }

        correlation_entry = {
            "sample": sample_name,
            "comparison": comp_name,
        }

        for method in corr_methods:
            if method not in metric_funcs:
                raise ValueError(f"Unknown correlation method: {method}")
            correlation_entry[f"{method}_corr"] = metric_funcs[method]()

        correlations.append(correlation_entry)

    return pd.DataFrame(correlations)


def plot_correlation(
    df: pd.DataFrame,
    corr_method: Literal["rmse", "pearson", "spearman", "mape", "ccc"] = "rmse",
) -> alt.Chart:
    """
    Create a heatmap with correlation values annotated as text.
    """
    num_replicates = df["sample"].nunique()
    num_comparisons = df["comparison"].nunique()
    chart_width = max(200, num_replicates * 100)
    chart_height = max(150, num_comparisons * 30)
    heatmap = (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X("sample:N", title="Sample"),
            y=alt.Y("comparison:N", title="Comparison"),
            color=alt.Color(
                f"{corr_method}_corr:Q",
                scale=alt.Scale(scheme="greens"),
                legend=alt.Legend(title=corr_method.capitalize()),
            ),
        )
        .properties(
            title=f"{corr_method.capitalize()} Correlation",
            width=chart_width,
            height=chart_height,
        )
    )

    text = (
        alt.Chart(df)
        .mark_text(baseline="middle", fontSize=12, color="white")
        .encode(
            x="sample:N",
            y="comparison:N",
            text=alt.Text(f"{corr_method}_corr:Q", format=".4f"),
        )
    )

    return heatmap + text


def normalize_sample_name(replicate_name: str) -> str:
    """
    Normalize replicate names to a sample key.
    - Illumina: remove trailing _REPxx
    - PacBio/Nanopore: normalize trailing replicate numbers
    """
    name = re.sub(r"_REP\d+$", "", replicate_name)  # Illumina
    name = re.sub(r"replicate\d+$", "replicate", name)  # PacBio/Nanopore
    return name or replicate_name


# ---------------- Main Workflow ---------------- #

replicate_dfs = {}
meth_caller_correlation_dfs = []
sample_correlation = []
corr_methods = snakemake.params["correlation_methods"]

################### COMPUTE CORRELATION ####################

for sample_file in snakemake.input["samples"]:
    # Get the name of the specific sample.
    replicate_name = sample_file.split("/")[-3]
    df = pd.read_parquet(sample_file, engine="pyarrow")

    meth_caller_correlation_dfs.append(
        compute_correlation(df, replicate_name, "meth_caller", corr_methods)
    )
    # Take the sample + replicate name and remove the specific replicate identifier
    samplename = normalize_sample_name(replicate_name)

    # Merge replicates if sample already exists (There are always two replicates, so the first one just gets inserted, the second one merged)
    if samplename in replicate_dfs:
        replicate_dfs[samplename] = pd.merge(
            replicate_dfs[samplename],
            df,
            on=["chromosome", "position"],
            how="inner",
            suffixes=("_rep1", "_rep2"),
        )
    else:
        replicate_dfs[samplename] = df
    if samplename.startswith("TrueMethylBS_HG002_LAB01"):
        print(replicate_name, file=sys.stderr)
        print(replicate_dfs[samplename].head(), file=sys.stderr)
        print(df.head(), file=sys.stderr)
# Save intermediate tables
with pd.HDFStore(snakemake.output["table"]) as store:
    for key, df in replicate_dfs.items():
        store[key] = df



# Calculate correlations across replicates
for sample_name, df in replicate_dfs.items():
    sample_correlation.append(
        compute_correlation(df, sample_name, "replicate", corr_methods)
    )

# Combine results
correlation_meth_callers = pd.concat(meth_caller_correlation_dfs, ignore_index=True)
correlation_replicates = pd.concat(sample_correlation, ignore_index=True)
################### PLOTTING ####################

charts_meth_callers_list = [
    plot_correlation(correlation_meth_callers, method) for method in corr_methods
]
charts_replicates_list = [
    plot_correlation(correlation_replicates, method) for method in corr_methods
]

charts_meth_callers = alt.hconcat(*charts_meth_callers_list).properties(
    title="Correlation of Varlo with Reference Tools"
)
charts_replicates = alt.hconcat(*charts_replicates_list).properties(
    title="Correlation of Replicates"
)

full_chart = alt.vconcat(charts_meth_callers, charts_replicates)

# Save final chart
full_chart.save(
    snakemake.output["plot"],
    embed_options={"actions": False},
    inline=False,
)
