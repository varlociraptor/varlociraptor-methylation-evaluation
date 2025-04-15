import pandas as pd
import altair as alt

# Redirect standard error to snakemake log file
sys.stderr = open(snakemake.log[0], "w")


pd.set_option("display.max_columns", None)

csv_files = snakemake.input
dfs = [pd.read_csv(file) for file in csv_files]
df_combined = pd.concat(dfs, ignore_index=True)
df_combined = df_combined.sort_values(by=["tool", "coverage", "recall"])

line_chart = (
    alt.Chart(df_combined[df_combined["tool"] == "varlo"])
    .mark_line()
    .encode(
        x="recall:Q",
        y=alt.Y(
            "precision:Q",
            scale=alt.Scale(
                domain=[
                    max(0, min(df_combined["precision"]) - 0.02),
                    min(1, max(df_combined["precision"]) + 0.02),
                ]
            ),
        ),
        color="tool:N",
        order=alt.Order("prob_pres_threshold", sort="ascending"),
    )
)

point_chart = (
    alt.Chart(df_combined)
    .mark_point()
    .encode(
        x="recall:Q",
        y="precision:Q",
        color="tool:N",
        tooltip=[
            "tool:N",
            "coverage:N",
            "number_sites:Q",
            "recall:Q",
            "precision:Q",
            "prob_pres_threshold:N",
        ],
    )
)

chart = (line_chart + point_chart).properties(
    width=800,
    height=400,
    title=f"Precision-Recall Curve for Coverage {df_combined['coverage'].iloc[0]}",
)

chart.save(snakemake.output[0], scale_factor=2.0)
