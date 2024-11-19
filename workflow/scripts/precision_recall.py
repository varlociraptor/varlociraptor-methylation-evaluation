import pandas as pd
import altair as alt

# Alle CSV-Dateien einlesen
csv_files = snakemake.input  # Beispiel-Pfad, je nach Struktur anpassen

# Alle CSV-Dateien in eine Liste von DataFrames einlesen
dfs = [pd.read_csv(file) for file in csv_files]

# Alle DataFrames zu einem großen DataFrame zusammenführen
df_combined = pd.concat(dfs, ignore_index=True)

# Daten sortieren nach "tool" und "coverage" für eine glattere Linie im Plot
df_combined = df_combined.sort_values(by=["tool", "coverage", "recall"])


print(df_combined)

# Altair Plot erstellen
chart = (
    alt.Chart(df_combined)
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
        color="tool:N",  # Verschiedene Farben für jedes Tool
        shape="coverage:N",  # Verschiedene Linienarten für unterschiedliche Coverages
        tooltip=[
            "tool:N",
            "coverage:N",
            "number_sites:Q",
            "recall:Q",
            "precision:Q",
        ],  # Zeige Details bei Hover an
    )
    .properties(
        width=800,
        height=400,
        title="Precision-Recall Curve for Different Tools and Coverages",
    )
)

# Speichern als HTML für interaktive Visualisierung mit Tooltips und Zoom
chart.save(snakemake.output[0], format=snakemake.params["plot_type"])
