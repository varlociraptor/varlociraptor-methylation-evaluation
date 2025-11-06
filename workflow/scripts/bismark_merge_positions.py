import pandas as pd
import pysam

# Read bedgraph
bedgraph = pd.read_csv(
    snakemake.input["bedgraph"],
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "perc", "meth", "unmeth"],
)
bedgraph["chrom"] = bedgraph["chrom"].astype(str)

# Read BCF
bcf = pysam.VariantFile(snakemake.input["candidates"][0])

results = []

for idx, row in bedgraph.iterrows():
    chrom = row["chrom"]
    pos_start = int(row["start"])
    pos_end = int(row["end"])

    records = list(bcf.fetch(chrom, pos_start, pos_end))
    if len(records) == 0:
        print(f"No candidate in BCF for {chrom}:{pos_start}-{pos_end}")
        continue
    bcf_position = records[0].pos
    print(records, chrom, pos_start, pos_end)
    total_meth = row["meth"]
    total_unmeth = row["unmeth"]
    # merged_percentage = (total_meth / (total_meth + total_unmeth)) * 100

    results.append(
        {
            "chrom": chrom,
            "pos": bcf_position,
            # "meth_pct_merged": round(merged_percentage, 2),
            # "coverage": int(total_meth + total_unmeth),
            "meth_counts": int(total_meth),
            "unmeth_counts": int(total_unmeth),
        }
    )

df = pd.DataFrame(results)

# Nach chrom und pos gruppieren, Summen von meth_counts und unmeth_counts berechnen
df_merged = df.groupby(["chrom", "pos"], as_index=False).agg(
    {"meth_counts": "sum", "unmeth_counts": "sum"}
)

# Neue Spalten berechnen
df_merged["coverage"] = df_merged["meth_counts"] + df_merged["unmeth_counts"]
df_merged["meth_pct_merged"] = round(
    df_merged["meth_counts"] / df_merged["coverage"] * 100, 2
)

# Als TSV speichern
df_merged.to_csv(snakemake.output[0], sep="\t", header=False, index=False)
