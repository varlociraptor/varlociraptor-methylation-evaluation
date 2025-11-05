import pandas as pd
import pysam

# Read bedgraph
bedgraph = pd.read_csv(
    snakemake.input["bedgraph"],
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "perc", "meth", "unmeth"],
)

# Read BCF
bcf_file = pysam.VariantFile(snakemake.input["candidates"][0])

results = []

for record in bcf_file:
    chrom = record.chrom
    pos = record.pos  # 1-based position
    # Find matching bedgraph positions
    mask = (bedgraph["chrom"] == chrom) & (
        (bedgraph["start"] == pos) | (bedgraph["end"] == pos)
    )
    matching = bedgraph[mask]

    if matching.empty:
        continue

    total_meth = matching["meth"].sum()
    total_unmeth = matching["unmeth"].sum()
    merged_percentage = (total_meth / (total_meth + total_unmeth)) * 100

    results.append(
        {
            "chrom": chrom,
            "pos_start": pos - 1,
            "pos_end": pos + 1,
            "meth_pct_merged": round(merged_percentage, 2),
            "coverage": int(total_meth + total_unmeth),
            "meth_counts": int(total_meth),
            "unmeth_counts": int(total_unmeth),
        }
    )

# Save merged results
merged_df = pd.DataFrame(results)
merged_df.to_csv(snakemake.output[0], sep="\t", header=False, index=False)
