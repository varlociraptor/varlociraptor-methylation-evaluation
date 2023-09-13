# Erweitern Sie die Wildcards und erstellen Sie eine Liste aller passenden Dateipfade
print(snakemake.input[0])
bedGraph_files = [snakemake.input[i] for i in range(len(snakemake.input))]

# expand("Pfad/zum/Verzeichnis/{wildcard}.bedGraph", wildcard=your_wildcard_pattern)
print(bedGraph_files)
# Öffnen Sie alle BedGraph-Dateien gleichzeitig
with [open(file_path, "r") for file_path in bedGraph_files] as bedGraph_files:
    # Hier können Sie den Code zum Einlesen und Verarbeiten der BedGraph-Dateien schreiben
    # bedGraph_files ist eine Liste von geöffneten Dateien

    # Zum Beispiel:
    for file in bedGraph_files:
        for line in file:
            # Verarbeiten Sie jede Zeile in der BedGraph-Datei
            pass
