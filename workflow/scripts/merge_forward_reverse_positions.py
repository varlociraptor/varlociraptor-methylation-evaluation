import pandas as pd
import pysam

# Read bedgraph
bedgraph = pd.read_csv(
    snakemake.input["bedgraph"],
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "perc", "col5", "col6"],
)

# Ensure chromosome is string
bedgraph["chrom"] = bedgraph["chrom"].astype(str)

# Drop rows with missing start or end
bedgraph = bedgraph.dropna(subset=["start", "end", "perc", "col5"])

# Bismark output has meth/unmeth column but no coverage column, BisSNP has only a coverage column
if not bedgraph["col6"].isna().all():
    # Bismark
    bedgraph["meth"] = bedgraph["col5"].astype(int)
    bedgraph["unmeth"] = bedgraph["col6"].astype(int)
else:
    # BisSNP
    bedgraph["coverage"] = bedgraph["col5"]
    bedgraph["meth"] = ((bedgraph["perc"] / 100) * bedgraph["coverage"]).round().astype(int)
    bedgraph["unmeth"] = (bedgraph["coverage"] - bedgraph["meth"]).astype(int)

# Read BCF
bcf = pysam.VariantFile(snakemake.input["candidates"][0])

results = []

for idx, row in bedgraph.iterrows():
    chrom = row["chrom"]
    pos_start = int(row["start"])
    pos_end = int(row["end"])

    # Fetch variants in BCF
    records = list(bcf.fetch(chrom, pos_start, pos_end))
    if len(records) == 0:
        continue  # no candidate, skip

    bcf_position = records[0].pos
    total_meth = row["meth"]
    total_unmeth = row["unmeth"]

    results.append({
        "chrom": chrom,
        "pos": bcf_position,
        "meth_counts": int(total_meth),
        "unmeth_counts": int(total_unmeth)
    })

# Create DataFrame
df = pd.DataFrame(results)

# Merge by chrom + pos
df_merged = df.groupby(["chrom", "pos"], as_index=False).agg({
    "meth_counts": "sum",
    "unmeth_counts": "sum"
})

# Compute coverage and merged methylation percentage
df_merged["coverage"] = df_merged["meth_counts"] + df_merged["unmeth_counts"]
df_merged["meth_pct_merged"] = round(
    df_merged["meth_counts"] / df_merged["coverage"] * 100, 2
)

# Save output
df_merged.to_csv(snakemake.output[0], sep="\t", header=False, index=False)
