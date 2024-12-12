import altair as alt
import pandas as pd
import numpy as np
import re
import os


def compute_bias(prob, format_values):
    probs = re.split("[=;]", prob)
    probs_values = [probs[i] for i in [1, 3, 5]]
    try:
        prob_values = [float(value) for value in probs_values]
    except:
        return "bias"
    min_value = min(prob_values)
    min_index = prob_values.index(min_value)

    bias_labels = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]

    if any(value != "." for value in format_values[6:12]):
        bias = bias_labels[format_values[6:12].index(".")]
    else:
        bias = "normal"

    if min_index != 0:
        bias = "normal"
    return bias


def compute_rmse(df, tool_col):
    squared_errors = (df[tool_col] - df["true_methylation"]) ** 2
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
        (df["varlo_methylation"] - df["true_methylation"]).abs().round().astype(int)
    )

    df["ref_distance"] = (
        (df["ref_tool_methylation"] - df["true_methylation"]).abs().round().astype(int)
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
            y="true_methylation",
            # color="bias:N",
            opacity=alt.value(0.3),
            tooltip=["chromosome:N", "position:N"],
        )
        .interactive()
    )
    scatter_ref = (
        alt.Chart(df)
        .mark_circle(color="red")
        .encode(
            x="ref_tool_methylation",
            y="true_methylation",
            # color="bias:N",
            opacity=alt.value(0.3),
            tooltip=["chromosome:N", "position:N"],
        )
        .interactive()
    )
    chart1 = (
        (scatter_meth + scatter_ref + line)
        .properties(
            width=400,
            height=400,
            title=alt.Title(
                f"Varlo and {ref_tool} vs. TrueMeth",
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
                f"Varlo and {ref_tool} vs. TrueMeth",
                subtitle=f"Varlo methylation (blue) in front of {ref_tool} methylation (red)",
            ),
        )
        .interactive()
    )

    # Beide Charts nebeneinander
    combined_chart = chart1 | chart2

    # Speichern
    combined_chart.save(output, scale_factor=2.0)


def debug_distances(df, output):
    print("Output: ", output)
    print(df)
    with open(output, "w") as out_file:
        df = df.sort_values(by="distance_tools", ascending=False)
        out_file.write(
            f"big_dist: {list(zip(df['chromosome'], df['position']))[:20]}\n"
        )

        df = df.sort_values(by="distance_tools", ascending=True)
        out_file.write(
            f"small_dist: {list(zip(df['chromosome'], df['position']))[:20]}\n"
        )


def read_truth_file(truth_file_path):
    # Lies TrueMeth-Daten
    df = []

    with open(truth_file_path, "r") as truth_file:
        for line in truth_file:
            if line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            chrom, start, end, meth_rate = (
                parts[0].replace("chr", ""),
                int(parts[1]),
                int(parts[2]),
                float(parts[3]),
            )
            position = (start + end) // 2
            df.append([chrom, position, meth_rate])

    return pd.DataFrame(df, columns=["chromosome", "position", "true_methylation"])


def read_tool_file(tool_file_path, file_name):
    # Allgemeine Struktur für die Datei
    df = []

    with open(tool_file_path, "r") as tool_file:

        for line in tool_file:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")

            # Spezifische Struktur für jedes Tool
            if file_name == "Varlociraptor":  # Varlociraptor
                chrom, pos, alternative, info_field, format_field, values = (
                    parts[0],
                    int(parts[1]),
                    parts[4],
                    parts[7],
                    parts[8],
                    parts[9].split(":"),
                )
                if alternative == "<METH>":
                    format_fields = format_field.split(":")
                    dp_index = format_fields.index("DP")
                    af_index = format_fields.index("AF")

                    meth_rate = float(values[af_index]) * 100
                    coverage = int(values[dp_index])

                    bias = compute_bias(info_field, values)

                    df.append([chrom, pos, meth_rate, coverage, bias])

            elif file_name == "methylDackel":
                chrom, start, end, meth_rate, coverage = (
                    parts[0],
                    int(parts[1]),
                    int(parts[2]),
                    float(parts[3]),
                    int(parts[4]) + int(parts[5]),
                )
                position = (start + end) // 2
                df.append([chrom, position, meth_rate, coverage, "normal"])

            elif file_name == "bsMap":
                chrom, position, meth_rate, coverage = (
                    parts[0],
                    int(parts[1]),
                    float(parts[4]) * 100,
                    int(parts[5]),
                )
                df.append([chrom, position, meth_rate, coverage, "normal"])

            elif file_name == "bismark":
                chrom, start, end, meth_rate = (
                    parts[0],
                    int(parts[1]),
                    int(parts[2]),
                    float(parts[3]),
                )
                position = (start + end) // 2
                df.append([chrom, position, meth_rate, None, "normal"])

            elif file_name == "bisSNP":
                chrom, position, meth_rate = (
                    parts[0],
                    int(parts[2]),
                    float(parts[3]),
                )
                df.append([chrom, position, meth_rate, None, "normal"])

            elif file_name == "modkit":
                details = parts[9].split()
                chrom, position, coverage, meth_rate = (
                    parts[0].removeprefix("chr"),
                    int(parts[2]),
                    int(details[0]),
                    float(details[1]),
                )
                df.append([chrom, position, meth_rate, None, "normal"])

            elif file_name == "pb_CpG_tools":
                chrom, position, meth_rate, coverage = (
                    parts[0],
                    int(parts[2]),
                    float(parts[3]),
                    int(parts[5]),
                )

                df.append([chrom, position, meth_rate, coverage, "normal"])

    # Erstelle DataFrame
    columns = [
        "chromosome",
        "position",
        "tool_methylation",
        "coverage",
        "bias",
    ]
    return pd.DataFrame(df, columns=columns)


pd.set_option("display.max_columns", None)
varlo_file = snakemake.input["tool"]
ref_file = snakemake.input["ref_tool"]
ref_tool = os.path.splitext(os.path.basename(ref_file))[0]


varlo_df = read_tool_file(varlo_file, "Varlociraptor")
ref_df = read_tool_file(ref_file, ref_tool)
truth_df = read_truth_file(snakemake.input["true_meth"][0])


df_temp = pd.merge(
    varlo_df[varlo_df["bias"] == "normal"],  # Gefiltertes varlo_df
    ref_df,  # ref_df
    on=["chromosome", "position"],  # Gemeinsame Spalten
    how="inner",  # Nur gemeinsame Einträge
)


df = pd.merge(
    df_temp,  # Ergebnis des ersten Merges
    truth_df,  # truth_df
    on=["chromosome", "position"],  # Gemeinsame Spalten
    how="inner",  # Nur gemeinsame Einträge
)

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

debug_distances(df, snakemake.output["distances"])
