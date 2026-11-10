import sys

import altair as alt
import polars as pl

METH_CALLER_TO_NAME = {
    "varlo_1.0": "Varlociraptor α = 1.0",
    "varlo_0.5": "Varlociraptor α = 0.5",
    "varlo_0.1": "Varlociraptor α = 0.1",
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

TOOL_COLORS = {
    "Bismark": "#D81B60",
    "BSMAPz": "#1E88E5",
    "BisSNP": "#FFC107",
    "MethylDackel": "#f0700e",
    "Modkit": "#B42CEA",
    "Pb-CpG-tools": "#8D9279",
    "Varlociraptor": "#028d77",
}

METRIC_TO_NAME = {
    "mape": "Relative discordance",
    "mae": "Absolute discordance",
    # "binary_discordance": " Dₛ",
}

# Points based on fewer loci are too noisy and are not shown
MIN_LOCI = 1000

pl.Config.set_tbl_cols(20)
pl.Config.set_tbl_rows(20)
alt.data_transformers.enable("vegafusion")

sys.stderr = open(snakemake.log[0], "w")


def read_coverage(coverage_file: str, rep_name: str) -> pl.DataFrame:
    return (
        pl.read_csv(
            coverage_file,
            separator="\t",
            has_header=False,
            new_columns=["chromosome", "pos_start", "pos_end", "coverage"],
        )
        .with_columns(
            ((pl.col("pos_end") + pl.col("pos_start")) / 2)
            .cast(pl.Int64)
            .alias("position")
        )
        .select(["chromosome", "position", "coverage"])
        .rename({"coverage": f"coverage_{rep_name}"})
    )


def load_methylation_data(meth_data_path: str, meth_callers: list[str]) -> pl.DataFrame:
    methylation_columns = [
        f"{caller}_methylation_{rep}"
        for caller in meth_callers
        for rep in ["rep1", "rep2"]
    ]
    return (
        pl.scan_parquet(meth_data_path)
        .select(["sample", "chromosome", "position", *methylation_columns])
        .with_columns(pl.col("chromosome").cast(pl.Int64))
        .collect()
    )


def join_coverage(
    df: pl.DataFrame, coverage_path_rep1: str, coverage_path_rep2: str
) -> pl.DataFrame:
    """
    Join the coverage of both replicates and compute their minimum per site.
    Sites without a coverage entry in either replicate get a null min coverage.
    """
    coverage_rep1 = read_coverage(coverage_path_rep1, "rep1")
    coverage_rep2 = read_coverage(coverage_path_rep2, "rep2")
    return (
        df.join(coverage_rep1, on=["chromosome", "position"], how="left")
        .join(coverage_rep2, on=["chromosome", "position"], how="left")
        .with_columns(
            pl.when(
                pl.col("coverage_rep1").is_not_null()
                & pl.col("coverage_rep2").is_not_null()
            )
            .then(pl.min_horizontal("coverage_rep1", "coverage_rep2"))
            .alias("min_coverage")
        )
    )


def discordance(
    pairs: pl.DataFrame, caller: str, alpha: float, min_coverage: int | None
) -> pl.DataFrame:
    """
    One row with the number of loci and the discordance between the replicates of a caller.
    Same formulas as compute_replicate_counts in prepare_plot_df.py.
    """
    rep1, rep2 = pl.col("rep1"), pl.col("rep2")
    abs_diff = (rep1 - rep2).abs()
    return pairs.select(
        N=pl.len(),
        mape=pl.when((rep1 == 0) & (rep2 == 0))
        .then(0.0)
        .otherwise(abs_diff / pl.max_horizontal(rep1, rep2))
        .mean(),
        mae=abs_diff.mean(),
        # binary_discordance=((rep1 > 0) != (rep2 > 0)).cast(pl.Float64).mean(),
        caller=pl.lit(caller),
        alpha=pl.lit(alpha),
        min_coverage=pl.lit(min_coverage, dtype=pl.Int64),
    )


def caller_info(df: pl.DataFrame, caller: str) -> pl.DataFrame:
    return df.select(
        pl.col(f"{caller}_methylation_rep1").alias("rep1"),
        pl.col(f"{caller}_methylation_rep2").alias("rep2"),
        "min_coverage",
    ).drop_nulls(["rep1", "rep2"])


def coverage_threshold(caller_df: pl.DataFrame, n_target: int) -> int:
    """
    Integer min coverage threshold t whose number of retained loci N(t) is closest to n_target.
    Ties go to the smaller t. Loci without coverage only pass t = 0.
    """
    # N(t) = number of loci with min_coverage >= t, i.e. floor(min_coverage) >= t
    loci_per_threshold = (
        caller_df.group_by(
            pl.col("min_coverage").fill_null(0).floor().cast(pl.Int64)
        )
        .len()
        .sort("min_coverage", descending=True)
        .with_columns(N=pl.col("len").cum_sum())
    )
    return loci_per_threshold.sort(
        (pl.col("N") - n_target).abs(), "min_coverage"
    ).item(0, "min_coverage")


def compute_discordance_vs_loci(
    df: pl.DataFrame, meth_callers: list[str]
) -> pl.DataFrame:
    """
    Varlociraptor is only filtered by its filter level α. Every α retains N_α loci. Every other
    caller gets the min coverage threshold whose number of retained loci is closest to N_α.
    """
    varlo_callers = [c for c in meth_callers if c.startswith("varlo_")]
    ref_callers = [c for c in meth_callers if not c.startswith("varlo_")]

    varlo_df = pl.concat(
        discordance(
            caller_info(df, caller),
            caller="varlociraptor",
            alpha=float(caller.removeprefix("varlo_")),
            min_coverage=None,
        )
        for caller in varlo_callers
    )

    frames = [varlo_df]
    for caller in ref_callers:
        caller_df = caller_info(df, caller)
        for alpha, n_alpha in varlo_df.select("alpha", "N").rows():
            threshold = coverage_threshold(caller_df, n_alpha)
            filtered = caller_df.filter(
                pl.col("min_coverage").fill_null(0) >= threshold
            )
            frames.append(
                discordance(filtered, caller, alpha, min_coverage=threshold)
            )
    return pl.concat(frames, how="diagonal_relaxed")


def plot_discordance_vs_loci(df: pl.DataFrame):
    df = df.filter(pl.col("N") > MIN_LOCI).with_columns(
        # Rightmost point of each curve (largest N) is drawn bigger
        pl.when(pl.col("N") == pl.col("N").max().over("sample", "metric", "tool_name"))
        .then(60)
        .otherwise(12)
        .alias("point_size")
    )
    color_domain = [
        tool for tool in TOOL_COLORS if tool in df["tool_name"].unique().to_list()
    ]
    color_range = [TOOL_COLORS[tool] for tool in color_domain]

    base = alt.Chart().encode(
        x=alt.X("N:Q", title="Retained CpG loci (N)"),
        y=alt.Y("value:Q", scale=alt.Scale(zero=False)),
        color=alt.Color(
            "tool_name:N",
            title="Caller",
            scale=alt.Scale(domain=color_domain, range=color_range),
        ),
    )
    line = base.mark_line(strokeWidth=1).encode(order=alt.Order("N:Q"))
    points = base.mark_point(filled=True, opacity=1).encode(
        size=alt.Size("point_size:Q", scale=None, legend=None)
    )

    return (
        alt.layer(line, points, data=df)
        .properties(width=200, height=200)
        .facet(
            row=alt.Row(
                "metric:N", title=None, sort=list(METRIC_TO_NAME.values())
            ),
            column=alt.Column("sample:N", title=None),
        )
        .resolve_scale(y="independent")
        .configure_header(labelFontSize=14, labelFontWeight="bold")
    )


meth_callers = snakemake.params["meth_callers"]
samples = snakemake.params["samples"]

platform = snakemake.params["platform"]


meth_df = load_methylation_data(snakemake.input["meth_data"], meth_callers)

sample_frames = []
for sample, coverage_path_rep1, coverage_path_rep2 in zip(
    samples, snakemake.input["coverages_rep1"], snakemake.input["coverages_rep2"]
):
    df = join_coverage(
        meth_df.filter(pl.col("sample") == sample), coverage_path_rep1, coverage_path_rep2
    )
    sample_frames.append(
        compute_discordance_vs_loci(df, meth_callers).with_columns(
            # PacBio and Nanopore only have the sample "dummy", so we show the platform instead
            pl.lit(
                platform if sample == "dummy" else sample.replace("_HG002_", "_")
            ).alias("sample")
        )
    )

df_discordance = (
    pl.concat(sample_frames)
    # Make a row per metric instead of a column per metric, so that we can facet the plot by metric. The new column is called "metric" and the values are in a column called "value".
    .unpivot(
        index=["sample", "caller", "alpha", "min_coverage", "N"],
        on=list(METRIC_TO_NAME.keys()),
        variable_name="metric",
    )
    # Rename metrics and callers to more readable names
    .with_columns(
        pl.col("metric").replace(METRIC_TO_NAME),
        pl.col("caller").replace(METH_CALLER_TO_NAME).alias("tool_name"),
    )
)


chart = plot_discordance_vs_loci(df_discordance)
# vl-convert rejects scale != 1 for PDF (vega_to_pdf path, taken when vegafusion pre-compiles the spec)
chart.save(
    snakemake.output[0],
    scale_factor=1 if snakemake.output[0].endswith(".pdf") else 2,
)
df_discordance.select(
    "sample", "caller", "tool_name", "alpha", "min_coverage", "N", "metric", "value"
).write_parquet(snakemake.output[1])
