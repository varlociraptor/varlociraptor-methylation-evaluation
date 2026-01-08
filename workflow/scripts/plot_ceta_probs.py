import pandas as pd
from pathlib import Path
from functools import reduce
import altair as alt
import polars as pl
import sys

# sys.stderr = open(snakemake.log[0], "w")

pl.Config.set_tbl_rows(-1)  # show all rows
pl.Config.set_tbl_cols(-1)  # show all columns

files = snakemake.input

name_to_method = {
    "MethylSeq_HG002_LAB01_REP01_with_prior": "MethylSeq with Prior",
    "EMSeq_HG002_LAB01_REP01_with_prior": "EMSeq with Prior",
    "MethylSeq_HG002_LAB01_REP01_no_prior": "MethylSeq no Prior",
    "EMSeq_HG002_LAB01_REP01_no_prior": "EMSeq no Prior",
    "ceta_multi": "Common",
    "ceta_multi_all": "Common all",
    "ceta_multi_both": "Common both",
    "ceta_multi_emseq": "Common emseq",
    "ceta_multi_untreated": "Common untreated",
    "untreated_with_prior": "Untreated with Prior",
    "untreated_no_prior": "Untreated no Prior",
}


colorblind_safe_palette = [
    "#034D40",
    "#078A72",
    "#0EC5A4",
    "#14F8CE",
    "#75F8E0",
    "#D81B60",
    "#E05387",
    "#386791",
    "#1E88E5",
    "#B89B46",
    "#FFC107",
]
dfs = []

for f in files:
    method = name_to_method.get(Path(f).parts[-3], "Unknown")

    # Load parquet
    df = pl.scan_parquet(f)
    df = df.with_columns(pl.lit(method).alias("method"))
    dfs.append(df)

df = pl.concat(dfs)

bin_size = 50
intervals = [i / bin_size for i in range(bin_size)]
print(intervals)
interval_labels = ["missing"] + [
    f"{i / bin_size} - { (i + 1) / bin_size}" for i in range(bin_size)
]

df = df.with_columns(
    pl.col("prob_present")
    .cut(intervals, labels=interval_labels, left_closed=True)
    .alias("prob_present_bin"),
    pl.col("prob_absent")
    .cut(intervals, labels=interval_labels, left_closed=True)
    .alias("prob_absent_bin"),
    pl.col("prob_artifact")
    .cut(intervals, labels=interval_labels, left_closed=True)
    .alias("prob_artifact_bin"),
)
print(df.collect().head())
# ---- Neue Filterspalte basierend auf cg_pos und c_to_t ----
# df = df.with_columns(
#     pl.when(pl.col("c_to_t") == True)
#     .then("c_to_t")
#     .otherwise("complete")
#     .alias("filter_option")
# )


def make_plot(df, category):
    df = df.collect().to_pandas()
    col = f"prob_{category}_bin"
    df[col] = df[col].cat.add_categories(["missing"]).fillna("missing")

    chart = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X(
                f"prob_{category}_bin:N",
                title=None,
                scale=alt.Scale(domain=interval_labels),
                axis=alt.Axis(labelAngle=-45),
            ),
            xOffset=alt.XOffset("method:N"),
            y=alt.Y("count()", title="Count"),
            color=alt.Color(
                "method:N",
                title="Method",
                scale=alt.Scale(range=colorblind_safe_palette),
                legend=alt.Legend(labelLimit=0),
            ),
            opacity=alt.condition(method_select, alt.value(1), alt.value(0.1)),
            tooltip=["method:N", f"prob_{category}_bin:N", "count()"],
        )
        .transform_filter(
            (filter_param == "complete")
            | ((filter_param == "cg") & (alt.datum.cg_pos))
            | ((filter_param == "c2t") & (alt.datum.c_to_t))
            | ((filter_param == "cg & c2t") & (alt.datum.cg_pos) & (alt.datum.c_to_t))
            | ((filter_param == "cg & !c2t") & (alt.datum.cg_pos) & (~alt.datum.c_to_t))
            | ((filter_param == "!cg & c2t") & (~alt.datum.cg_pos) & (alt.datum.c_to_t))
            | (
                (filter_param == "!cg & !c2t")
                & (~alt.datum.cg_pos)
                & (~alt.datum.c_to_t)
            )
        )
        .transform_filter(
            alt.FieldOneOfPredicate(field="method", oneOf=name_to_method.values())
        )
        .add_params(filter_param, method_select)
        .add_params(method_select)
        .properties(title=f"Distribution of prob_{category}")
    )

    return chart


# method_radio = alt.binding_radio(options=list(name_to_method.values()), name="Method")
# method_select = alt.selection_point(fields=["method"], bind=method_radio)

method_select = alt.selection_point(fields=["method"])

filter_param = alt.param(
    name="Filter",
    bind=alt.binding_select(
        options=[
            "complete",
            "cg",
            "c2t",
            "cg & c2t",
            "cg & !c2t",
            "!cg & c2t",
            "!cg & !c2t",
        ],
        name="Filter: ",
    ),
    value="complete",
)

plot_present = make_plot(df, "present")
plot_absent = make_plot(df, "absent")
plot_artifact = make_plot(df, "artifact")

heatmap_plots = alt.vconcat(plot_present, plot_absent, plot_artifact).resolve_scale(
    y="shared"
)

heatmap_plots.save(snakemake.output[0], embed_options={"actions": False}, inline=False)
