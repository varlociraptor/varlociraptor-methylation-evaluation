import sys
import re
import numpy as np
import pandas as pd
import altair as alt
from typing import Literal

# Redirect stderr to Snakemake log file
sys.stderr = open(snakemake.log[0], "w")

# Pandas display options (useful for debugging)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)


def compute_ccc(x: pd.Series, y: pd.Series) -> float:
    """
    Compute the Concordance Correlation Coefficient (CCC) between two variables.

    Reference:
        Lin, L. I.-K. (1989). "A concordance correlation coefficient to evaluate reproducibility."
        Biometrics, 45(1), 255–268.
    """
    mean_x, mean_y = x.mean(), y.mean()
    var_x, var_y = x.var(ddof=0), y.var(ddof=0)
    pearson_r = x.corr(y, method="pearson")

    numerator = 2 * pearson_r * np.sqrt(var_x * var_y)
    denominator = var_x + var_y + (mean_x - mean_y) ** 2
    return numerator / denominator


def compute_correlation(
    df: pd.DataFrame,
    sample_name: str,
    mode: Literal["meth_caller", "replicate"] = "meth_caller",
    corr_methods: list[str] = ["rmse"],
) -> pd.DataFrame:
    """
    Compute correlations and error metrics for methylation data.

    Parameters:
        df (pd.DataFrame): Data containing methylation values.
        sample_name (str): Sample or replicate identifier.
        mode (str):
            - 'meth_caller': Compare varlociraptor against other tools.
            - 'replicate': Compare replicate 1 vs replicate 2 for the same tool.
        corr_methods (list[str]): List of correlation or distance metrics to compute.

    Returns:
        pd.DataFrame: Summary of computed correlation metrics.
    """
    results = []

    for meth_caller in snakemake.params["meth_callers"]:
        if mode == "meth_caller":
            # Skip Varlociraptor as it's the baseline reference
            if meth_caller == "varlo":
                continue
            col1, col2 = "varlo_methylation", f"{meth_caller}_methylation"
            comparison = f"varlo vs. {meth_caller}"

        elif mode == "replicate":
            col1, col2 = (
                f"{meth_caller}_methylation_rep1",
                f"{meth_caller}_methylation_rep2",
            )
            if col1 not in df.columns or col2 not in df.columns:
                continue
            comparison = meth_caller

        else:
            raise ValueError(f"Invalid mode: {mode}")

        # Drop NA rows for fair comparison
        clean_df = df.dropna(subset=[col1, col2])
        if clean_df.empty:
            continue

        # Define available metrics
        metric_funcs = {
            "pearson": lambda: clean_df[col1].corr(clean_df[col2], method="pearson"),
            "spearman": lambda: clean_df[col1].corr(clean_df[col2], method="spearman"),
            "rmse": lambda: np.sqrt(((clean_df[col1] - clean_df[col2]) ** 2).mean())
            / 100.0,
            "mape": lambda: (
                np.where(
                    (clean_df[col1] == 0) & (clean_df[col2] == 0),
                    0,
                    np.abs(clean_df[col1] - clean_df[col2])
                    / clean_df[[col1, col2]].max(axis=1),
                ).mean()
                * 100.0
            ),
            "ccc": lambda: compute_ccc(clean_df[col1], clean_df[col2]),
        }

        entry = {"sample": sample_name, "comparison": comparison}

        for method in corr_methods:
            if method not in metric_funcs:
                raise ValueError(f"Unknown correlation method: {method}")
            entry[f"{method}_corr"] = metric_funcs[method]()

        results.append(entry)

    return pd.DataFrame(results)


def normalize_sample_name(replicate_name: str) -> str:
    """
    Normalize replicate names to a sample identifier.
    """
    name = re.sub(r"_REP\d+$", "", replicate_name)  # Illumina pattern
    name = re.sub(
        r"minimal\d+$", "minimal", name
    )  # PacBio/ Nanopore / multi-sample pattern
    name = re.sub(r"minimal\d+$", "minimal", name)  # minimal toy example
    return name or replicate_name


replicate_dfs = {}
correlation_dfs = []
corr_methods = snakemake.params["correlation_methods"]

# ---- Compute correlations across all input samples ---- #

for sample_file in snakemake.input["samples"]:
    # Extract replicate name from file path
    replicate_name = sample_file.split("/")[-3]
    df = pd.read_parquet(sample_file, engine="pyarrow")

    # Compute correlation across tools for this replicate
    correlation_dfs.append(
        compute_correlation(
            df, replicate_name, mode="meth_caller", corr_methods=corr_methods
        )
    )

    # Normalize sample name (combine replicates)
    sample_name = normalize_sample_name(replicate_name)

    if sample_name in replicate_dfs:
        # Merge replicate 1 and replicate 2 data
        replicate_dfs[sample_name] = pd.merge(
            replicate_dfs[sample_name],
            df,
            on=["chromosome", "position"],
            how="inner",
            suffixes=("_rep1", "_rep2"),
        )
    else:
        replicate_dfs[sample_name] = df

# ---- Save intermediate tables ---- #

with pd.HDFStore(snakemake.output["table"]) as store:
    for key, data in replicate_dfs.items():
        store[key] = data
