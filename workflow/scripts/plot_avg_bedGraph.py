import pickle
import numpy as np
import altair as alt
import pandas as pd
import random

candidates = snakemake.input.candidates
bedGraph_files = snakemake.input.bedgraphs

# bedGraph_files = [snakemake.input[i] for i in range(len(snakemake.input["bedgraphs"]))]

positions_to_graphs = {}
bedGraph_entry = {}


all_keys = set()


with open(candidates, "r") as file:
    for line in file:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        chrom, pos = "chr" + parts[0], parts[1]
        positions_to_graphs[(chrom, int(pos))] = []


for i, file_path in enumerate(bedGraph_files):
    file_keys = set()
    print(i, file_path)
    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            chrom, start, end, methylation, meth_reads, unmeth_reads = parts
            pos = int((int(start) + int(end)) / 2)
            coverage = meth_reads + unmeth_reads

            current_data = positions_to_graphs.get((chrom, pos), None)
            if current_data is not None:
                positions_to_graphs[(chrom, pos)].append(
                    (i, coverage, methylation))
            else:
                print(f"Key ({chrom}, {pos}) not found, data not appended.")

data = []
for (chrom, pos), values in positions_to_graphs.items():
    for i, coverage, methylation in values:
        data.append({'chrom': chrom, 'pos': pos, 'i': i,
                    'coverage': coverage, 'methylation': float(methylation)})
df = pd.DataFrame(data)
unique_positions = df['pos'].unique()
random_positions = random.sample(
    list(unique_positions), min(20, len(unique_positions)))
df = df[df['pos'].isin(random_positions)]

# Sortieren Sie den Datensatz nach dem durchschnittlichen Methylierungswert pro Position

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
print(df)
# Sortieren Sie den DataFrame nach der durchschnittlichen Methylierung pro Position
df['mean_methylation'] = df.groupby('pos')['methylation'].transform('mean')
df = df.sort_values(by=['mean_methylation', "pos"])
print(df)

# 1. Graph: Coverage
chart1 = alt.Chart(df).mark_circle().encode(
    x=alt.X('pos:O', title='Positions', axis=alt.Axis(labelAngle=0),
            sort=alt.EncodingSortField(field="mean_methylation", order='ascending')),
    y=alt.Y('coverage:Q', title='Coverage'),
    color=alt.Color('i:N', scale=alt.Scale(scheme='category20'), legend=None)
).properties(width=400, title='Coverage per position')

# 2. Graph: Methylierung
chart2 = alt.Chart(df).mark_circle().encode(
    x=alt.X('pos:O', title='Positions', axis=alt.Axis(labelAngle=0),
            sort=alt.EncodingSortField(field="mean_methylation", order='ascending')),
    y=alt.Y('methylation:Q', title='Methylation'),
    color=alt.Color('i:N', scale=alt.Scale(scheme='category20'), legend=None)
).properties(width=400, title='Methylation per position')

# Grafiken anzeigen
chart1.save(snakemake.output["cov"], scale_factor=2.0)
chart2.save(snakemake.output["meth"], scale_factor=2.0)
