import altair as alt
import pandas as pd
import numpy as np
import re
import os
from scipy.stats import linregress


# Redirect standard error to snakemake log file
sys.stderr = open(snakemake.log[0], "w")


def compute_rmse(df, tool_col):
    # df_rmse = df[df[tool_col] > 0]
    df_rmse = df
    squared_errors = (df_rmse[tool_col] - df_rmse["true_methylation_x"]) ** 2
    mean_squared_error = squared_errors.mean()
    return np.sqrt(mean_squared_error)


def linear_regression_summary(x, y):
    result = linregress(x, y)
    slope = round(result.slope, 2)
    intercept = round(result.intercept, 2)
    r_value = round(result.rvalue, 2)
    r_squared = round(result.rvalue**2, 2)
    return slope, intercept, r_value, r_squared


def scatter_plot(df, ref_tool, output):
    rmse_varlo = round(compute_rmse(df, "varlo_methylation"), 2)
    rmse_ref = round(compute_rmse(df, "ref_tool_methylation"), 2)
    try:
        m, b, r, r2 = linear_regression_summary(
            df["ref_tool_methylation"], df["varlo_methylation"]
        )
    except:
        m = "Not valid"
        b = "Not valid"
        r = "Not valid"
        r2 = "Not valid"

    line = (
        alt.Chart(pd.DataFrame({"x": [0, 100], "y": [0, 100]}))
        .mark_line(color="green")
        .encode(x="x:Q", y="y:Q")
    )

    scatter_varlo = (
        alt.Chart(df)
        .mark_circle(color="blue")
        .encode(
            x=alt.X("varlo_methylation", title="Varlo Methylation"),
            y=alt.Y("true_methylation_x", title=f"True Methylation"),
            # color="bias:N",
            opacity=alt.value(0.3),
            tooltip=["chromosome:N", "position:N", "coverage_x:N"],
        )
        .interactive()
    )
    scatter_ref = (
        alt.Chart(df)
        .mark_circle(color="red")
        .encode(
            x=alt.X("ref_tool_methylation", title=f"{ref_tool} Methylation"),
            y=alt.Y("true_methylation_x", title=f"True Methylation"),
            # color="bias:N",
            opacity=alt.value(0.3),
            tooltip=["chromosome:N", "position:N", "coverage_y:N"],
        )
        .interactive()
    )
    scatter_both = (
        alt.Chart(df)
        .mark_circle(color="green")
        .encode(
            x=alt.X("ref_tool_methylation", title=f"{ref_tool} Methylation"),
            y=alt.Y("varlo_methylation", title="Varlo Methylation"),
            # color="bias:N",
            opacity=alt.value(0.3),
            tooltip=["chromosome:N", "position:N", "coverage_y:N"],
        )
        .interactive()
    )
    regression_line = (
        alt.Chart(df)
        .transform_regression(
            "ref_tool_methylation", "varlo_methylation", method="linear"
        )
        .mark_line(color="black", strokeDash=[5, 5])
        .encode(x="ref_tool_methylation:Q", y="varlo_methylation:Q")
    )

    chart1 = (
        (scatter_varlo + scatter_ref + line)
        .properties(
            width=400,
            height=400,
            title=alt.Title(
                f"Varlo (RMSE {rmse_varlo}) and {ref_tool} (RMSE {rmse_ref}) vs. TrueMeth",
                subtitle=f"{ref_tool} methylation (red) in front of Varlo methylation (blue)",
            ),
        )
        .interactive()
    )
    chart2 = (
        (scatter_ref + scatter_varlo + line)
        .properties(
            width=400,
            height=400,
            title=alt.Title(
                f"Varlo (RMSE {rmse_varlo}) and {ref_tool} (RMSE {rmse_ref}) vs. TrueMeth",
                subtitle=f"Varlo methylation (blue) in front of {ref_tool} methylation (red)",
            ),
        )
        .interactive()
    )
    chart3 = (
        (scatter_both + line + regression_line)
        .properties(
            width=400,
            height=400,
            title=alt.Title(
                f"Varlo (RMSE {rmse_varlo}) vs. {ref_tool} (RMSE {rmse_ref})",
                subtitle=f"R = {r}, R² = {r2}, y = {m}·x + {b}",
            ),
        )
        .interactive()
    )

    # Beide Charts nebeneinander
    combined_chart = chart1 | chart2 | chart3

    # Speichern
    combined_chart.save(output, scale_factor=2.0)


pd.set_option("display.max_columns", None)
varlo_file = snakemake.input["varlo"]
ref_file = snakemake.input["ref_tool"]
ref_tool = os.path.splitext(os.path.basename(ref_file))[0].capitalize()

varlo_df = pd.read_parquet(varlo_file, engine="pyarrow")
# varlo_df = varlo_df[
#     varlo_df["prob_present"] >= float(snakemake.params["prob_pres_threshhold"])
# ]
ref_df = pd.read_parquet(ref_file, engine="pyarrow")


df = pd.merge(
    varlo_df[varlo_df["bias"] == "normal"],  # Gefiltertes varlo_df
    ref_df,  # ref_df
    on=["chromosome", "position"],  # Gemeinsame Spalten
    how="inner",  # Nur gemeinsame Einträge
)


#######################################
#  Outer Merge: Alle Zeilen aus beiden DataFrames
merged_outer = pd.merge(
    varlo_df,
    ref_df,
    on=["chromosome", "position"],
    how="outer",
    indicator=True,  # Fügt eine Spalte "_merge" hinzu
)


# Zeilen filtern, die nur in einem der beiden DataFrames enthalten sind
not_in_merged = merged_outer[merged_outer["_merge"] != "both"]

#####################################

df.rename(
    columns={
        "tool_methylation_x": "varlo_methylation",
        "tool_methylation_y": "ref_tool_methylation",
    },
    inplace=True,
)

scatter_plot(df, ref_tool, snakemake.output["scatter_plot"][0])
