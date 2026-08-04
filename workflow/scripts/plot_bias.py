import sys

import altair as alt
import polars as pl

# Logging
sys.stderr = open(snakemake.log[0], "w")

pl.Config.set_tbl_rows(100)
pl.Config.set_tbl_cols(200)
alt.data_transformers.enable("vegafusion")


BIAS_COLS = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]
INFO_COLS = ["DP", "AF"]
KEEP_COLS = INFO_COLS + BIAS_COLS
VARLO_COLS = INFO_COLS + ["SAOBS", "SROBS", "OBS", "OOBS"] + BIAS_COLS + ["AFD"]

BIAS_LABELS = {
    "SB": "Strand Bias",
    "ROB": "Read Orientation Bias",
    "RPB": "Read Position Bias",
    "SCB": "Soft-clipped Bias",
    "HE": "Haplotype Error",
    "ALB": "Alt Locus Bias",
}

REPLICATE_LABELS = {
    "DP_rep1": "replicate 1",
    "DP_rep2": "replicate 2",
}


def split_varlo_format(df: pl.DataFrame, rep: str, fdr: str) -> pl.DataFrame:
    """Split the colon-separated Varlociraptor FORMAT field into named columns."""
    col = f"varlo_{fdr}_format_{rep}"
    n_fields = len(VARLO_COLS)

    fields = df.select(pl.col(col).str.splitn(":", n_fields).alias("f")).unnest("f")
    fields.columns = [f"{name}_{rep}" for name in VARLO_COLS[: fields.width]]

    keep = [f"{name}_{rep}" for name in KEEP_COLS if f"{name}_{rep}" in fields.columns]
    return fields.select(keep)


def classify_bias_expr() -> pl.Expr:
    """Assign each variant to a bias category based on both replicates."""
    r1, r2 = pl.col("rep1_has_bias"), pl.col("rep2_has_bias")
    af1, af2 = pl.col("AF_rep1"), pl.col("AF_rep2")

    return (
        pl.when(r1 & r2)
        .then(pl.lit("Bias both reps"))
        .when((r1 & (af2 == 0)) | (r2 & (af1 == 0)))
        .then(pl.lit("Bias, AF = 0"))
        .when((r1 & (af2 > 0)) | (r2 & (af1 > 0)))
        .then(pl.lit("Bias, AF > 0"))
        .otherwise(pl.lit("No bias"))
    )


def build_bias_dataframe(df: pl.DataFrame, fdr: str) -> pl.DataFrame:
    """Build a long-format bias-analysis dataframe from Varlociraptor data."""
    df_r1 = split_varlo_format(df, "rep1", fdr)
    df_r2 = split_varlo_format(df, "rep2", fdr)
    base = pl.concat([df.select(["chromosome", "position", "sample"]), df_r1, df_r2], how="horizontal")
    print(base)
    bias_fields = [
        f"{bias}_{rep}" for bias in BIAS_COLS for rep in ("rep1", "rep2") if f"{bias}_{rep}" in base.columns
    ]
    if not bias_fields:
        return pl.DataFrame(schema={"chromosome": pl.Utf8, "position": pl.Int64})

    base = base.filter(pl.all_horizontal([pl.col(c).is_not_null() for c in bias_fields]))
    base = base.filter(pl.any_horizontal([pl.col(c) != "." for c in bias_fields]))
    if base.is_empty():
        return pl.DataFrame()

    base = base.with_columns(
        pl.col("AF_rep1", "AF_rep2").cast(pl.Float64),
        pl.col("DP_rep1", "DP_rep2").cast(pl.Int64),
    )
    base = base.with_columns(
        pl.any_horizontal([pl.col(f"{bias}_rep1") != "." for bias in BIAS_COLS]).alias("rep1_has_bias"),
        pl.any_horizontal([pl.col(f"{bias}_rep2") != "." for bias in BIAS_COLS]).alias("rep2_has_bias"),
    )
    base = base.with_columns(classify_bias_expr().alias("category"))
    print(base)
    id_vars = ["chromosome", "position", "sample", "AF_rep1", "AF_rep2", "DP_rep1", "DP_rep2", "category"]
    long = base.unpivot(index=id_vars, on=bias_fields, variable_name="bias_var", value_name="bias_value")
    print(long)
    long = long.filter(pl.col("bias_value") != ".")

    long = long.with_columns(
        pl.col("bias_var").str.replace(r"_rep[12]$", "").alias("bias_type"),
        pl.col("bias_var").str.extract(r"_(rep[12])$", 1).alias("replicate"),
    )
    long = long.with_columns(pl.col("bias_type").replace(BIAS_LABELS).alias("bias_type_label"))
    print(long)
    return long


def bias_plots(df_long: pl.DataFrame, fdr: str, platform_label: str) -> alt.HConcatChart:
    """Create bias, allele-frequency, and depth plots from long-format data."""
    bias_chart = (
        alt.Chart(df_long)
        .mark_bar()
        .encode(
            x=alt.X(
                "category:N",
                axis=alt.Axis(labelAngle=-45),
                title=None,
                scale=alt.Scale(domain=["Bias both reps", "Bias, AF = 0", "Bias, AF > 0"]),
            ),
            y="count():Q",
            color=alt.Color(
                "bias_type_label:N",
                scale=alt.Scale(
                    domain=df_long["bias_type_label"].unique().to_list(),
                    range=["#D81B60", "#1E88E5"],
                ),
                title="Bias Type",
            ),
            tooltip=["category", "count()", "bias_type_label"],
        )
        .properties(title=platform_label)
    )

    df_af = (
        df_long.filter(pl.col("category") == "Bias, AF > 0")
        .with_columns(pl.max_horizontal("AF_rep1", "AF_rep2").round(2).alias("AF"))
        # .filter((pl.col("DP_rep1") <= 500) & (pl.col("DP_rep2") <= 500))
        .with_columns(
            pl.when(pl.col("AF_rep1") == 0).then(pl.col("DP_rep1")).otherwise(pl.col("DP_rep2")).alias("DP_bias"),
            pl.when(pl.col("AF_rep1") > 0).then(pl.col("DP_rep1")).otherwise(pl.col("DP_rep2")).alias("DP_AF"),
        )
    )
    print(df_af.filter((pl.col("DP_bias") > 1000) | (pl.col("DP_AF") > 1000)))
    dp_scatter = (
        alt.Chart(df_af)
        .mark_circle()
        .encode(
            x=alt.X("DP_bias:Q", title="Depth associated with bias"),
            y=alt.Y("DP_AF:Q", title="Depth associated with AF > 0"),
            tooltip=["DP_bias:Q", "DP_AF:Q"],
        )
        .properties(title="Depth scatter plot for 'Bias, AF > 0'")
    )

    af_chart = (
        alt.Chart(df_af)
        .mark_bar(color="#05AA8F")
        .encode(
            x=alt.X("AF:Q", bin=alt.Bin(step=0.05), title="Allele Frequency"),
            y="count():Q",
            tooltip=["AF:Q", "count():Q"],
        )
        .properties(title="AF at one sided biased loci")
    )

    df_dp = df_af.unpivot(
        index=["chromosome", "position"],
        on=["DP_rep1", "DP_rep2"],
        variable_name="replicate",
        value_name="DP",
    ).with_columns(pl.col("replicate").replace(REPLICATE_LABELS))

    dp_chart = (
        alt.Chart(df_dp)
        .mark_bar()
        .encode(
            x=alt.X("DP:Q", bin=alt.Bin(maxbins=50), title="Depth"),
            y="count():Q",
            color=alt.Color(
                "replicate:N",
                scale=alt.Scale(domain=list(REPLICATE_LABELS.values()), range=["#FFC107", "#004D40"]),
            ),
            tooltip=["DP:Q", "count():Q"],
        )
        .properties(title="Coverage distributions")
    )


    return (
        alt.hconcat(bias_chart, af_chart, dp_chart, dp_scatter)
        .resolve_scale(color="independent")
        .properties(title=f"FDR {fdr}")
    )


def empty_plot(fdr: str) -> alt.Chart:
    """Placeholder chart shown when there is no bias data for a given FDR."""
    return (
        alt.Chart(pl.DataFrame({"msg": [f"No bias data for FDR {fdr}"]}))
        .mark_text(size=20)
        .encode(text="msg:N")
    )


samples = snakemake.params["sample"]
if isinstance(samples, str):
    samples = [samples]

df = pl.read_parquet(snakemake.input[0])
df = df.filter(pl.col("sample").is_in(samples))

platform = snakemake.params["platform"]
platform_label = "Illumina" if platform == "Illumina_pe" else platform

all_charts = []
df_long = pl.DataFrame()

for fdr in snakemake.params["fdrs"]:
    cols = [
        "chromosome",
        "position",
        f"varlo_{fdr}_format_rep1",
        f"varlo_{fdr}_format_rep2",
        f"varlo_{fdr}_methylation_rep1",
        f"varlo_{fdr}_methylation_rep2",
        "sample"
    ]
    df_subset = df.select(cols)
    print(df_subset)
    df_long = build_bias_dataframe(df_subset, fdr)
    chart = empty_plot(fdr) if df_long.is_empty() else bias_plots(df_long, fdr, platform_label)
    all_charts.append(chart)

final_chart = alt.vconcat(*all_charts)
final_chart.save(snakemake.output[0])
df_long.write_parquet(snakemake.output[1])
