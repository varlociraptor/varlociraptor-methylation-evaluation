import os
import pandas as pd


def read_truth_file(truth_file_path, pos_big_dist, pos_small_dist):
    # Lies TrueMeth-Daten

    with open(truth_file_path, "r") as truth_file:
        lines_big_dist = {}
        lines_small_dist = {}
        for line in truth_file:
            if line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            chrom, start, end = (
                parts[0].replace("chr", ""),
                int(parts[1]),
                int(parts[2]),
            )
            position = (start + end) // 2
            
            if (chrom, position) in pos_big_dist:
                lines_big_dist[(chrom, position)] = line
            if (chrom, position) in pos_small_dist:
                lines_small_dist[(chrom, position)] = line


    return lines_big_dist, lines_small_dist



def read_tool_file(tool_file_path, file_name, pos_big_dist, pos_small_dist):
    lines_big_dist = {}
    lines_small_dist = {}

    with open(tool_file_path, "r") as varlo_file:
        for line in varlo_file:
            chrom, position = None, None
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")

            if file_name == "Varlociraptor":
                chrom, position = (
                    str(parts[0]),
                    int(parts[1]),
                )

            elif file_name == "methylDackel":
                chrom, start, end = (
                    str(parts[0]),
                    int(parts[1]),
                    int(parts[2]),
                )
                position = (start + end) // 2

            elif file_name == "bsMap":
                chrom, position = (
                    str(parts[0]),
                    int(parts[1]),
                )

            elif file_name == "bismark":
                chrom, start, end = (
                    str(parts[0]),
                    int(parts[1]),
                    int(parts[2]),
                )

            elif file_name == "bisSNP":
                chrom, position = (
                    str(parts[0]),
                    int(parts[2]),
                )
            elif file_name == "modkit":
                chrom, position = (
                    str(parts[0]).removeprefix("chr"),
                    int(parts[2]),
                )

            elif file_name == "pb_CpG_tools":
                chrom, position = (
                    str(parts[0]),
                    int(parts[2]),
                )

            if (chrom, position) in pos_big_dist:
                lines_big_dist[(chrom, position)] = line
            if (chrom, position) in pos_small_dist:
                lines_small_dist[(chrom, position)] = line

            # Weitere File-Typen können hier ähnlich hinzugefügt werden, wie bereits im Original-Skript

    return lines_big_dist, lines_small_dist


def read_distances_file(distances_file_path):
    pos_big_dist = []
    pos_small_dist = []

    with open(distances_file_path, "r") as dist_file:

        for line in dist_file:
            line = line.strip()
            if line.startswith("big_dist"):
                list_str = line.split(":")[1].strip()
                big_positions = eval(list_str)
                pos_big_dist = big_positions
            elif line.startswith("small_dist"):
                list_str = line.split(":")[1].strip()
                small_positions = eval(list_str)
                pos_small_dist = small_positions

    return pos_big_dist, pos_small_dist





def extract_lines(varlo_file, ref_tool_file, truth_file, distances_file, output_file):
    # Extrahiere Positionen aus der distances-Datei
    pos_big_dist, pos_small_dist = read_distances_file(distances_file)

    # Bestimme den Dateinamen basierend auf der Eingabedatei
    ref_tool_name = os.path.splitext(os.path.basename(ref_tool_file))[0]

    # Liest die Tool-Datei und extrahiert Zeilen mit den passenden Positionen
    lines_big_dist_varlo, lines_small_dist_varlo = read_tool_file(
        varlo_file, "Varlociraptor", pos_big_dist, pos_small_dist
    )

    lines_big_dist_tool, lines_small_dist_tool = read_tool_file(
        ref_tool_file, ref_tool_name, pos_big_dist, pos_small_dist
    )
    
    lines_big_dist_truth, lines_small_dist_truth = read_truth_file(
        truth_file, pos_big_dist, pos_small_dist
    )
    # Speichere die extrahierten Zeilen in der Ausgabedatei
    with open(output_file, "w") as out_file:
        out_file.write("big_dist_lines:\n")
        for pos in pos_big_dist:
            out_file.write(lines_big_dist_varlo.get(pos, "Not found"))
            out_file.write(lines_big_dist_tool.get(pos, "Not found"))
            out_file.write(lines_big_dist_truth.get(pos, "Not found"))
            out_file.write("\n")

        out_file.write("\n\n\n\n\nsmall_dist_lines:\n")
        for pos in pos_small_dist:
            out_file.write(lines_small_dist_varlo.get(pos, "Not found"))
            out_file.write(lines_small_dist_tool.get(pos, "Not found"))
            out_file.write(lines_small_dist_truth.get(pos, "Not found"))
            
            out_file.write("\n")


# Hauptfunktion für die Verarbeitung
if __name__ == "__main__":
    varlo_file = snakemake.input["tool"]
    ref_tool_file = snakemake.input["ref_tool"]
    truth_file = snakemake.input["true_meth"][0]
    distances_file = snakemake.input["distances"]
    output_file = snakemake.output["distances"]

    extract_lines(varlo_file, ref_tool_file, truth_file, distances_file, output_file)
