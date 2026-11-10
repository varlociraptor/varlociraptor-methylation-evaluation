import sys

import altair as alt
import numpy as np
import pandas as pd

# Logging
sys.stderr = open(snakemake.log[0], "w")

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
alt.data_transformers.enable("vegafusion")


BIAS_COLS = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]
INFO_COLS = ["DP", "AF"]
KEEP_COLS = INFO_COLS + BIAS_COLS
# FORMAT fields of the Varlociraptor calls (branch fdr-filter-bins). The parquet only stores the values.
VARLO_FORMAT = ["DP", "AF", "SAOBS", "SROBS", "OBS", "OOBS", "SB", "ROB", "RPB", "SCB", "HE", "ALB", "FEB", "AFD"]

BIAS_LABELS = {
    "SB": "Strand Bias",
    "ROB": "Read Orientation Bias",
    "RPB": "Read Position Bias",
    "SCB": "Soft‑clipped Bias",
    "HE": "Haplotype Error",
    "ALB": "Alt Locus Bias",
}

BIAS_COLORS = {
    "Strand Bias": "#D81B60",
    "Alt Locus Bias": "#1E88E5",
    "Read Orientation Bias": "#FFC107",
    "Read Position Bias": "#004D40",
    "Soft‑clipped Bias": "#8E44AD",
    "Haplotype Error": "#FF7F0E",
}

REPLICATE_LABELS = {
    "DP_rep1": "replicate 1",
    "DP_rep2": "replicate 2",
}


def split_varlo_format(df: pd.DataFrame, rep: str, fdr: str) -> pd.DataFrame:
    """Extract colon-separated Varlociraptor FORMAT fields."""
    col = f"varlo_{fdr}_format_{rep}"
    fields = df[col].str.split(":", expand=True)
    if fields.shape[1] != len(VARLO_FORMAT):
        raise ValueError(f"Expected {len(VARLO_FORMAT)} FORMAT fields in {col}, got {fields.shape[1]}")
    fields.columns = [f"{c}_{rep}" for c in VARLO_FORMAT]
    return fields[[f"{c}_{rep}" for c in KEEP_COLS]]


def classify_bias(df):
    """Vectorized bias category assignment."""
    r1, r2 = df["rep1_has_bias"], df["rep2_has_bias"]
    af1, af2 = df["AF_rep1"], df["AF_rep2"]

    return np.select(
        [
            r1 & r2,
            (r1 & (af2 == 0)) | (r2 & (af1 == 0)),
            (r1 & (af2 > 0)) | (r2 & (af1 > 0)),
        ],
        ["Bias both reps", "Bias, AF = 0", "Bias, AF > 0"],
        default="No bias",
    )


def build_bias_dataframe(df: pd.DataFrame, fdr: str) -> pd.DataFrame:
    """Build bias-analysis dataframe from Varlociraptor data with one row per locus (sample, chromosome, position)."""
    df_r1 = split_varlo_format(df, "rep1", fdr)
    df_r2 = split_varlo_format(df, "rep2", fdr)

    base = pd.concat(
        [df[["sample", "chromosome", "position"]].reset_index(drop=True), df_r1, df_r2],
        axis=1,
    )

    bias_fields = [f"{b}_{r}" for b in BIAS_COLS for r in ("rep1", "rep2")]
    base = base[base[bias_fields].notna().all(axis=1)]
    base = base[(base[bias_fields] != ".").any(axis=1)]

    base[["AF_rep1", "AF_rep2"]] = base[["AF_rep1", "AF_rep2"]].astype(float)
    base[["DP_rep1", "DP_rep2"]] = base[["DP_rep1", "DP_rep2"]].astype(int)

    base["rep1_has_bias"] = base[[f"{b}_rep1" for b in BIAS_COLS]].ne(".").any(axis=1)
    base["rep2_has_bias"] = base[[f"{b}_rep2" for b in BIAS_COLS]].ne(".").any(axis=1)

    base["category"] = classify_bias(base)
    if base.empty:
        return pd.DataFrame()
    long = base.melt(
        id_vars=[
            "sample",
            "chromosome",
            "position",
            "AF_rep1",
            "AF_rep2",
            "DP_rep1",
            "DP_rep2",
            "category",
        ],
        value_vars=bias_fields,
        var_name="bias_var",
        value_name="bias_value",
    )

    long = long[long["bias_value"] != "."]
    long[["bias_type", "replicate"]] = long["bias_var"].str.rsplit(
        pat="_", n=1, expand=True
    )

    long["bias_type_label"] = long["bias_type"].map(BIAS_LABELS)

    # Count loci, not bias entries: a locus with a bias in both replicates is one locus. Varlociraptor calls at most one
    # bias per replicate; if the replicates differ, the bias type of replicate 1 is used.
    locus = ["sample", "chromosome", "position"]
    labels = (
        long.sort_values("replicate", kind="stable")
        .groupby(locus)["bias_type_label"]
        .first()
    )
    loci = (
        long.drop_duplicates(locus)
        .drop(columns=["bias_var", "bias_value", "bias_type", "replicate", "bias_type_label"])
        .merge(labels.reset_index(), on=locus)
    )
    return loci


def bias_plots(df_long: pd.DataFrame, fdr: str, platform_label: str):
    """Create bias, AF, and DP plots from long-format data."""
    present_types = set(df_long["bias_type_label"])
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
            y=alt.Y("count():Q", title="Number of loci"),
            color=alt.Color(
                "bias_type_label:N",
                # Only the bias types present, so the legend has no unused entries
                scale=alt.Scale(
                    domain=[t for t in BIAS_COLORS if t in present_types],
                    range=[c for t, c in BIAS_COLORS.items() if t in present_types],
                ),
                title="Bias Type",
                # legend=None if platform_label != "Nanopore" else alt.Legend(),
            ),
            tooltip=["category", "count()", "bias_type_label"],
        )
    ).properties(width=150, height=200, title=f"VarLociraptor FDR={fdr}")

    # AF plot
    # df_af = df_long[df_long["category"] == "Bias, AF > 0"].assign(
    #     AF=lambda d: d[["AF_rep1", "AF_rep2"]].max(axis=1).round(2)
    # )
    # df_af = df_af[(df_af["DP_rep1"] <= 500) & (df_af["DP_rep2"] <= 500)]

    # af_chart = (
    #     alt.Chart(df_af)
    #     .mark_bar(color="#05AA8F")
    #     .encode(
    #         x=alt.X("AF:Q", bin=alt.Bin(step=0.05), title="Allele Frequency"),
    #         y="count():Q",
    #         tooltip=["AF:Q", "count():Q"],
    #     )
    #     .properties(title="AF at one sided biased loci")
    # )

    # # DP plot
    # df_dp = df_af.melt(
    #     id_vars=["chromosome", "position"],
    #     value_vars=["DP_rep1", "DP_rep2"],
    #     var_name="replicate",
    #     value_name="DP",
    # ).assign(replicate=lambda d: d["replicate"].map(REPLICATE_LABELS))
    # dp_chart = (
    #     alt.Chart(df_dp)
    #     .mark_bar()
    #     .encode(
    #         x=alt.X("DP:Q", bin=alt.Bin(maxbins=50), title="Depth"),
    #         y="count():Q",
    #         color=alt.Color(
    #             "replicate:N",
    #             scale=alt.Scale(
    #                 domain=list(REPLICATE_LABELS.values()),
    #                 range=["#FFC107", "#004D40"],
    #             ),
    #         ),
    #         tooltip=["DP:Q", "count():Q"],
    #     )
    #     .properties(title="Coverage distributions")
    # )

    # final_chart = (
    #     alt.hconcat(bias_chart, af_chart, dp_chart)
    #     .resolve_scale(color="independent")
    #     .properties(title=f"FDR {fdr}")
    # )
    return bias_chart


def empty_plot(fdr):
    """Create an empty plot with a message."""
    chart = (
        alt.Chart(pd.DataFrame({"msg": [f"No bias data for FDR {fdr}"]}))
        .mark_text(size=20)
        .encode(text="msg:N")
    )
    return chart


samples = snakemake.params["sample"]
if isinstance(samples, str):
    samples = [samples]
df = pd.read_parquet(snakemake.input[0], engine="pyarrow")
df = df[df["sample"].isin(samples)]

platform = snakemake.params["platform"]
platform_label = "Illumina" if platform == "Illumina_pe" else platform
all_charts = []
# Long-format bias data of every α, marked by column "fdr" (used by concat_bias.py)
bias_dfs = []
for fdr in snakemake.params["fdrs"]:
    cols = [
        "chromosome",
        "position",
        "sample",
        f"varlo_{fdr}_format_rep1",
        f"varlo_{fdr}_format_rep2",
        f"varlo_{fdr}_methylation_rep1",
        f"varlo_{fdr}_methylation_rep2",
    ]

    df_subset = df[cols]

    df_long = build_bias_dataframe(df_subset, fdr)
    if df_long.empty:
        chart = empty_plot(fdr)
    else:
        chart = bias_plots(df_long, fdr, platform_label)
        bias_dfs.append(df_long.assign(fdr=str(fdr)))
    all_charts.append(chart)
final_chart = alt.hconcat(*all_charts).properties(title=f"{platform_label}")
final_chart.save(snakemake.output[0])
df_all = (
    pd.concat(bias_dfs, ignore_index=True)
    if bias_dfs
    else pd.DataFrame(columns=["fdr"])
)
df_all.to_parquet(snakemake.output[1])
