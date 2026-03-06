import sys

import altair as alt
import numpy as np
import polars as pl

sys.stderr = open(snakemake.log[0], "w")
pl.Config.set_fmt_str_lengths(1000)
pl.Config.set_tbl_rows(100)
pl.Config.set_tbl_cols(100)
alt.data_transformers.enable("vegafusion")

# cast chromosome and position to string and int respectively
coverage = pl.read_csv(
    snakemake.input["coverage"],
    separator="\t",
    has_header=False,
    new_columns=["chromosome", "pos_start", "pos_end", "coverage"],
).with_columns(
    pl.col("chromosome").cast(pl.Utf8),
    pl.col("pos_start").cast(pl.Int64),
    pl.col("pos_end").cast(pl.Int64),
    pl.col("coverage").cast(pl.Int64),
    (pl.col("pos_start") + 1).alias("position"),
)

meth_data = pl.read_parquet(snakemake.input["meth_data"])

df = coverage.join(
    meth_data,
    on=["chromosome", "position"],
    how="outer",
).select(
    "chromosome",
    "position",
    "coverage",
    "varlo_0.01_methylation",
    "varlo_0.1_methylation",
    "varlo_1.0_methylation",
)
df = df.drop_nulls(subset=["coverage"])


# For each coverage in [0, 2, 5, 10, 20, 50] and varlo threshhold, compute the percentage of CpG sites that have at least that coverage and are retained in the final results (i.e. are present in the meth_data).
coverage_thresholds = [0, 1, 2, 3, 4, 5, 7, 10, 15, 20, 50, 100]
results = []
fdrs = [0.01, 0.1, 1.0]
for fdr in fdrs:
    meth_col = f"varlo_{fdr}_methylation"
    df_red = df.select("chromosome", "position", "coverage", pl.col(meth_col))

    df_red = df_red.drop_nulls()
    total_sites = df_red.shape[0]

    for cov in coverage_thresholds:
        retained_sites = df_red.filter((pl.col("coverage") >= cov))
        percentage_retained = (
            (retained_sites.shape[0] / total_sites * 100) if total_sites > 0 else np.nan
        )
        results.append(
            {
                "coverage_threshold": cov,
                "fdr": fdr,
                "percentage_retained": percentage_retained,
                "absolute_retained": retained_sites.shape[0],
            }
        )
# Plot the results as a line plot with coverage_threshold on the x-axis and percentage_retained on the y-axis, with a separate line for each fdr.
results_df = pl.DataFrame(results)


def plot_chart(y_axis="percentage_retained"):

    chart = (
        alt.Chart(results_df.to_pandas())
        .mark_line()
        .encode(
            x="coverage_threshold:O",
            y=f"{y_axis}:Q",
            color="fdr:N",
        )
        .properties(title=f"Coverage retained plot")
    )
    return chart


chart_relative = plot_chart("percentage_retained")
chart_absolute = plot_chart("absolute_retained")
chart = alt.hconcat(chart_relative, chart_absolute)
chart.save(snakemake.output[0])
