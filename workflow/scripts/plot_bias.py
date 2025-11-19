import pandas as pd
import altair as alt
import numpy as np
import sys

# Logging
sys.stderr = open(snakemake.log[0], "w")

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 1000)
alt.data_transformers.enable("vegafusion")


# --------------------------------------------------------------------
# CONSTANTS
# --------------------------------------------------------------------

BIAS_COLS = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]
INFO_COLS = ["DP", "AF"]
KEEP_COLS = INFO_COLS + BIAS_COLS
VARLO_COLS = INFO_COLS + ["SAOBS", "SROBS", "OBS", "OOBS"] + BIAS_COLS + ["AFD"]

BIAS_LABELS = {
    "SB": "Strand Bias",
    "ALB": "Alt Locus Bias",
}

REPLICATE_LABELS = {
    "DP_rep1": "replicate 1",
    "DP_rep2": "replicate 2",
}


def split_varlo_format(df: pd.DataFrame, rep: str, fdr: str) -> pd.DataFrame:
    """Extract colon-separated Varlociraptor FORMAT fields."""
    col = f"varlo_{fdr}_format_{rep}"
    fields = df[col].str.split(":", expand=True)
    fields.columns = [f"{c}_{rep}" for c in VARLO_COLS[: fields.shape[1]]]
    return fields[[f"{c}_{rep}" for c in KEEP_COLS if f"{c}_{rep}" in fields]]


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
    """Build long-format bias-analysis dataframe from Varlociraptor data."""
    df_r1 = split_varlo_format(df, "rep1", fdr)
    df_r2 = split_varlo_format(df, "rep2", fdr)

    base = pd.concat(
        [df[["chromosome", "position"]].reset_index(drop=True), df_r1, df_r2],
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

    long = base.melt(
        id_vars=[
            "chromosome",
            "position",
            "AF_rep1",
            "AF_rep2",
            "DP_rep1",
            "DP_rep2",
            "category",
        ],
        value_vars=[f"{b}_{r}" for b in BIAS_COLS for r in ("rep1", "rep2")],
        var_name="bias_var",
        value_name="bias_value",
    )

    long = long[long["bias_value"] != "."]

    long[["bias_type", "replicate"]] = long["bias_var"].str.rsplit(
        pat="_", n=1, expand=True
    )

    long["bias_type_label"] = long["bias_type"].map(BIAS_LABELS)

    return long


def make_plots(df_long: pd.DataFrame, fdr: str, platform_label: str):
    """Create bias, AF, and DP plots from long-format data."""
    # Bias category plot

    bias_chart = (
        alt.Chart(df_long)
        .mark_bar()
        .encode(
            x=alt.X("category:N", axis=alt.Axis(labelAngle=-45), title=None),
            y="count():Q",
            color=alt.Color(
                "bias_type_label:N",
                scale=alt.Scale(
                    domain=list(BIAS_LABELS.values()),
                    range=["#D81B60", "#1E88E5"],
                ),
                title="Bias Type",
            ),
            tooltip=["category", "count()", "bias_type_label"],
        )
    ).properties(title=f"{platform_label}")

    # AF plot
    df_af = df_long[df_long["category"] == "Bias, AF > 0"].assign(
        AF=lambda d: d[["AF_rep1", "AF_rep2"]].max(axis=1).round(2)
    )
    df_af = df_af[(df_af["DP_rep1"] <= 500) & (df_af["DP_rep2"] <= 500)]

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

    # DP plot
    df_dp = df_af.melt(
        id_vars=["chromosome", "position"],
        value_vars=["DP_rep1", "DP_rep2"],
        var_name="replicate",
        value_name="DP",
    ).assign(replicate=lambda d: d["replicate"].map(REPLICATE_LABELS))
    dp_chart = (
        alt.Chart(df_dp)
        .mark_bar()
        .encode(
            x=alt.X("DP:Q", bin=alt.Bin(maxbins=50), title="Depth"),
            y="count():Q",
            color=alt.Color(
                "replicate:N",
                scale=alt.Scale(
                    domain=list(REPLICATE_LABELS.values()),
                    range=["#FFC107", "#004D40"],
                ),
            ),
            tooltip=["DP:Q", "count():Q"],
        )
        .properties(title="Coverage distributions")
    )

    final_chart = (
        alt.hconcat(bias_chart, af_chart, dp_chart)
        .resolve_scale(color="independent")
        .properties(title=f"FDR {fdr}")
    )
    return final_chart


# --------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------

with pd.HDFStore(snakemake.input[0], mode="r", locking=False) as store:
    dfs = {key.strip("/"): store[key] for key in store.keys()}

samples = snakemake.params["sample"]
if isinstance(samples, str):
    samples = [samples]

sample_df = pd.concat([dfs[s] for s in samples], ignore_index=True)

platform = snakemake.params["platform"]
platform_label = "Illumina" if platform == "Illumina_pe" else platform

all_charts = []

for fdr in snakemake.params["fdrs"]:
    cols = [
        "chromosome",
        "position",
        f"varlo_{fdr}_format_rep1",
        f"varlo_{fdr}_format_rep2",
        f"varlo_{fdr}_methylation_rep1",
        f"varlo_{fdr}_methylation_rep2",
    ]

    df_subset = sample_df[cols]
    df_long = build_bias_dataframe(df_subset, fdr)
    chart = make_plots(df_long, fdr, platform_label)
    all_charts.append(chart)

final_chart = alt.vconcat(*all_charts)
final_chart.save(snakemake.output[0])
