import pandas as pd

# Bismark output a methylation rate for the forward and reverse strand separately without giving strand information. As a result we need to merge the positions manually. For this we compare the CpG position with our CpG candidate file and merge the positions accordingly.

bedgraph = pd.read_csv(
    snakemake.input["bedgraph"],
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "perc", "meth", "unmeth"],
)


vcf = pd.read_csv(
    snakemake.input["candidates"][0],
    comment="#",  
    sep="\t",
    header=None,
    names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
)


results = []

for _, row in vcf.iterrows():
    chrom = row["CHROM"]
    pos = row["POS"]
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

merged_df = pd.DataFrame(results).to_csv(
    snakemake.output[0], sep="\t", header=False, index=False
)
