from turtle import color, width

import altair as alt
import pandas as pd
import polars as pl

sys.stderr = open(snakemake.log[0], "w")

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)
pd.set_option("display.expand_frame_repr", False)
pd.set_option("display.max_rows", None)
pl.Config.set_tbl_cols(20)
pl.Config.set_tbl_rows(20)

sys.stderr = open(snakemake.log[0], "w")

alt.data_transformers.enable("vegafusion")

meth_caller_to_name = {
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

tool_colors = {
    "Bismark": "#D81B60",
    "BSMAPz": "#1E88E5",
    "BisSNP": "#FFC107",
    "MethylDackel": "#f0700e",
    "Modkit": "#B42CEA",
    "Pb-CpG-tools": "#8D9279",
    "Varlociraptor α = 0.1": "#004D40",
    "Varlociraptor α = 0.05": "#126e5f",
    "Varlociraptor α = 0.01": "#05AA8F",
}

df_illumina = pl.read_parquet(snakemake.input["illumina"]).with_columns(
    platform=pl.lit("Illumina")
)
df_pacbio = pl.read_parquet(snakemake.input["pacbio"]).with_columns(
    platform=pl.lit("PacBio")
)
df_nanopore = pl.read_parquet(snakemake.input["nanopore"]).with_columns(
    platform=pl.lit("Nanopore")
)
df = pl.concat([df_illumina, df_pacbio, df_nanopore])

df = df.filter(pl.col("min_coverage_bin") <= 80)

color_domain = sorted(df["tool_name"].unique())
color_range = [tool_colors[t] for t in color_domain]
line_plot_min_cov_vs_count = (
    alt.Chart(df)
    .mark_line()
    .encode(
        x=alt.X("min_coverage_bin:Q", title="Min Coverage"),
        y=alt.Y("cumulative_fraction:Q", title="Fraction of called CpG loci (N)"),
        color=alt.Color(
            "tool_name",
            title="Caller",
            scale=alt.Scale(
                domain=color_domain,
                range=color_range,
            ),
        ),
        strokeWidth=alt.value(1),
        column=alt.Column("platform", title=None),
    )
).properties(width=200, height=200).configure_header(labelFontSize=14, labelFontWeight="bold",)

line_plot_min_cov_vs_count.save(snakemake.output[0], scale_factor=2)
