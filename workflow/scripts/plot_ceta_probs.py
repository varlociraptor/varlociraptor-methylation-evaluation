import pandas as pd
from pathlib import Path
from functools import reduce
import altair as alt
import polars as pl
import sys

# sys.stderr = open(snakemake.log[0], "w")

pl.Config.set_tbl_rows(-1)  # show all rows
pl.Config.set_tbl_cols(-1)  # show all columns
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
files = snakemake.input

name_to_method = {
    "MethylSeq_HG002_LAB01_REP01_with_prior": "MethylSeq with Prior",
    "EMSeq_HG002_LAB01_REP01_with_prior": "EMSeq with Prior",
    # "MethylSeq_HG002_LAB01_REP01_no_prior": "MethylSeq no Prior",
    # "EMSeq_HG002_LAB01_REP01_no_prior": "EMSeq no Prior",
    # "ceta_multi": "Common",
    "ceta_multi_all": "Common all",
    "ceta_multi_emseq_untreated": "Common emseq untreated",
    "ceta_multi_emseq_methylseq": "Common emseq methylseq",
    # "ceta_multi_untr": "Common untreated",
    # "ceta_multi_no_untreated": "Common no untreated",
    "untreated_with_prior": "Untreated with Prior",
    # "untreated_no_prior": "Untreated no Prior",
}


colorblind_safe_palette = [
    "#034D40",
    "#078A72",
    "#0EC5A4",
    # "#14F8CE",
    # "#26FAEC",
    "#D81B60",
    # "#E05387",
    # "#386791",
    "#1E88E5",
    # "#B89B46",
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
bin_size = 20
intervals = [i / bin_size for i in range(bin_size)]
interval_labels = ["no coverage"] + [
    f"{i/ bin_size} - {(i+1)/ bin_size}" for i in range(bin_size)
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
# ---- Neue Filterspalte basierend auf cg_pos und c_to_t ----
# df = df.with_columns(
#     pl.when(pl.col("c_to_t") == True)
#     .then("c_to_t")
#     .otherwise("complete")
#     .alias("filter_option")
# )


def make_plot(df, category):
    col = f"prob_{category}_bin"
    print(df.filter(pl.col("method") == "MethylSeq with Prior").filter(pl.col("c_to_t") == True).head(200))
    grouped = (
        df
        .group_by([
            "method",
            col,
            "c_to_t",
            "cg_pos",
        ])
        .count()
        .collect()
    )

    # ---- Alle möglichen Kategorien erzeugen ----
    all_methods = pl.DataFrame({"method": list(name_to_method.values())})
    all_bins = pl.DataFrame({col: interval_labels}).with_columns(
        pl.col(col).cast(pl.Categorical)
    )
    all_c_to_t = pl.DataFrame({"c_to_t": [True, False]})
    all_cg_pos = pl.DataFrame({"cg_pos": [True, False]})

    # kartesisches Produkt
    full_grid = (
        all_methods
        .join(all_bins, how="cross")
        .join(all_c_to_t, how="cross")
        .join(all_cg_pos, how="cross")
    )


    # fehlende Kombinationen mit 0 auffüllen
    df = (
        full_grid
        .join(grouped, on=["method", col, "c_to_t", "cg_pos"], how="left")
        .with_columns(pl.col("count").fill_null(0))
        .to_pandas()
    )


    df[col] = df[col].fillna("no coverage")
    base = (
        alt.Chart(df)
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
            alt.FieldOneOfPredicate(field="method", oneOf=list(name_to_method.values()))
        )
    )
    print("Base chart created", category)

    base_line = base.transform_filter(
        alt.datum[col] != "no coverage"
    )

    line = base_line.mark_line(strokeWidth=2).encode(
        x=alt.X(
            f"{col}:N",
            title=None,
            scale=alt.Scale(domain=interval_labels),
            axis=alt.Axis(labelAngle=-45),
        ),
        y=alt.Y("sum(count):Q", title="Count"),
        color=alt.Color(
            "method:N",
            title="Method",
            scale=alt.Scale(range=colorblind_safe_palette),
            legend=alt.Legend(labelLimit=0),
        ),
        opacity=alt.condition(method_select, alt.value(1), alt.value(0.1)),
    )
    print("Line chart created", category)

    points = base_line.mark_point(size=15, filled=True).encode(
        x=alt.X(f"{col}:N", scale=alt.Scale(domain=interval_labels)),
        y=alt.Y("sum(count):Q"),
        color=alt.Color("method:N", scale=alt.Scale(range=colorblind_safe_palette)),
        opacity=alt.condition(method_select, alt.value(1), alt.value(0.1)),
        tooltip=["method:N", f"{col}:N", "sum(count):Q"],
    )

    points_no_cov = base.transform_filter(
        alt.datum[col] == "no coverage"
    ).mark_point(
        size=15,
        filled=True,
        shape="circle"
    ).encode(
        x=alt.X(f"{col}:N", scale=alt.Scale(domain=interval_labels)),
        y=alt.Y("sum(count):Q"),
        color=alt.Color(
            "method:N",
            scale=alt.Scale(range=colorblind_safe_palette),
        ),
        opacity=alt.condition(method_select, alt.value(1), alt.value(0.1)),
        tooltip=["method:N", f"{col}:N", "sum(count):Q"],
    )
    print("Points created", category)

    return (
        (line + points + points_no_cov)
        .add_params(filter_param, method_select)
        .properties(title=f"Distribution of prob_{category}")
    )


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

# df = df.collect().to_pandas()


plot_present = make_plot(df, "present")
plot_absent = make_plot(df, "absent")
plot_artifact = make_plot(df, "artifact")

heatmap_plots = alt.vconcat(plot_present, plot_absent, plot_artifact).resolve_scale(
    y="shared"
)
print("Charts created, saving...")

heatmap_plots.save(snakemake.output[0], embed_options={"actions": False})
