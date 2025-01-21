import altair as alt
import pandas as pd
import numpy as np
import re
import os


def compute_rmse(df, tool_col):
    squared_errors = (df[tool_col] - df["true_methylation_x"]) ** 2
    mean_squared_error = squared_errors.mean()
    return np.sqrt(mean_squared_error)


def compute_distances(df):
    """
    Plots the distribution of distances between predicted methylation values and truth methylation values for two tools.

    Args:
        tool1_name: Name of the first tool.
        tool2_name: Name of the second tool.
        true_meth: Name of the truth methylation column.
        tool1_values: List of methylation values predicted by the first tool.
        tool2_values: List of methylation values predicted by the second tool.
        true_values: List of truth methylation values.
        output: Path to the output file.
    """
    df["varlo_distance"] = (
        (df["varlo_methylation"] - df["true_methylation_x"]).abs().round().astype(int)
    )

    df["ref_distance"] = (
        (df["ref_tool_methylation"] - df["true_methylation_x"])
        .abs()
        .round()
        .astype(int)
    )

    df["distance_tools"] = (df["varlo_methylation"] - df["ref_tool_methylation"]).abs()

    melted_data = df.melt(
        value_vars=["ref_distance", "varlo_distance"],
        var_name="Tool",
        value_name="Distanz",
    )
    # Tool-Namen lesbarer machen
    melted_data["Tool"] = melted_data["Tool"].replace(
        {"ref_distance": "Reference", "varlo_distance": "Varlo"}
    )

    return df, melted_data

    # melted_data = melted_data[melted_data['Distanz'] <= 30]
    # Base chart with shared y-axis encoding


def distance_plot(melted_df, df, ref_tool, output):
    rmse_varlo = round(compute_rmse(df, "varlo_methylation"), 2)
    rmse_ref = round(compute_rmse(df, "ref_tool_methylation"), 2)
    base = (
        alt.Chart(melted_df)
        .mark_line()
        .encode(
            x=alt.X("Distanz:Q", title="Distance in %"),
            y=alt.Y("count():Q", title="Count"),
            color=alt.Color("Tool:N", legend=alt.Legend(title="Tool")),
        )
    )

    # Line chart
    line = base.mark_line().properties(
        title=f"Varlociraptor rmse: {str(rmse_varlo)} / {ref_tool} rmse: {str(rmse_ref)}",
    )

    # Point chart
    # points = base.mark_point()
    line.save(
        output,
    )


def scatter_plot(df, ref_tool, output):
    rmse_varlo = round(compute_rmse(df, "varlo_methylation"), 2)
    rmse_ref = round(compute_rmse(df, "ref_tool_methylation"), 2)

    line = (
        alt.Chart(pd.DataFrame({"x": [0, 100], "y": [0, 100]}))
        .mark_line(color="green")
        .encode(x="x:Q", y="y:Q")
    )

    scatter_meth = (
        alt.Chart(df)
        .mark_circle(color="blue")
        .encode(
            x="varlo_methylation",
            y="true_methylation_x",
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
            x="ref_tool_methylation",
            y="true_methylation_x",
            # color="bias:N",
            opacity=alt.value(0.3),
            tooltip=["chromosome:N", "position:N", "coverage_y:N"],
        )
        .interactive()
    )
    chart1 = (
        (scatter_meth + scatter_ref + line)
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
        (scatter_ref + scatter_meth + line)
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

    # Beide Charts nebeneinander
    combined_chart = chart1 | chart2

    # Speichern
    combined_chart.save(output, scale_factor=2.0)


pd.set_option("display.max_columns", None)
varlo_file = snakemake.input["varlo"]
ref_file = snakemake.input["ref_tool"]
ref_tool = os.path.splitext(os.path.basename(ref_file))[0]

varlo_df = pd.read_parquet(varlo_file, engine="pyarrow")
varlo_df = varlo_df[
    varlo_df["prob_present"] > float(snakemake.params["prob_pres_threshhold"])
]
ref_df = pd.read_parquet(ref_file, engine="pyarrow")

df = pd.merge(
    varlo_df[varlo_df["bias"] == "normal"],  # Gefiltertes varlo_df
    ref_df,  # ref_df
    on=["chromosome", "position"],  # Gemeinsame Spalten
    how="inner",  # Nur gemeinsame Einträge
)

print("Datensatz1: ", df)

df.rename(
    columns={
        "tool_methylation_x": "varlo_methylation",
        "tool_methylation_y": "ref_tool_methylation",
    },
    inplace=True,
)

scatter_plot(df, ref_tool, snakemake.output["scatter_plot"][0])

df, melted_data = compute_distances(df)
distance_plot(melted_data, df, ref_tool, snakemake.output["plot"][0])
