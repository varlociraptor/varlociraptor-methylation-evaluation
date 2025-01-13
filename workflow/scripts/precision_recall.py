import pandas as pd
import altair as alt

# CSV-Datei einlesen

# Alle CSV-Dateien einlesen
csv_files = snakemake.input  # Beispiel-Pfad, je nach Struktur anpassen

# Alle CSV-Dateien in eine Liste von DataFrames einlesen
dfs = [pd.read_csv(file) for file in csv_files]

# Alle DataFrames zu einem großen DataFrame zusammenführen
df_combined = pd.concat(dfs, ignore_index=True)

# Daten sortieren nach "tool" und "coverage" für eine glattere Linie im Plot
df_combined = df_combined.sort_values(by=["tool", "coverage", "recall"])

print(snakemake.output[0])

# Altair Plot erstellen: Linie und Punkte
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
        color="tool:N",  # Optional: Farbe der Punkte nach Tool
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

# Beide Charts kombinieren
chart = (line_chart + point_chart).properties(
    width=800,
    height=400,
    title=f"Precision-Recall Curve for Coverage {df_combined['coverage'].iloc[0]}",
)

# Speichern als HTML für interaktive Visualisierung
chart.save(snakemake.output[0], scale_factor=2.0)
