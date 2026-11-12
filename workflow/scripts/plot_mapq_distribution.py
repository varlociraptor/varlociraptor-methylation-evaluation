import sys

import altair as alt
import pandas as pd
import pysam

sys.stderr = open(snakemake.log[0], "w")
alt.data_transformers.enable("vegafusion")

chromosome = snakemake.params["chromosome"]

genome_label, chrom_label = "Genome", f"Chr{chromosome}"

# Number of primary, mapped reads per (region, MAPQ), in the whole alignment and on the chromosome
counts = {}
with pysam.AlignmentFile(snakemake.input[0], "rb") as bam:
    tid = bam.get_tid(chromosome)
    # until_eof: read the unsorted BAM sequentially, it has no index
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        mapq = read.mapping_quality
        counts[genome_label, mapq] = counts.get((genome_label, mapq), 0) + 1
        if read.reference_id == tid:
            counts[chrom_label, mapq] = counts.get((chrom_label, mapq), 0) + 1

df = pd.DataFrame(
    [(region, mapq, count) for (region, mapq), count in counts.items()],
    columns=["region", "mapq", "count"],
)
# Relative frequency per region; the number of reads n per region goes into the legend
total = df.groupby("region")["count"].transform("sum")
df["fraction"] = df["count"] / total
# df["region"] += " (n = " + total.map("{:,}".format) + ")"

legend_order = df["region"].unique().tolist()
chart = (
    alt.Chart(df)
    .mark_bar()
    .encode(
        x=alt.X("mapq:O", title="MAPQ", axis=alt.Axis(values=list(range(0, 61, 10)), labelAngle=0)),
        xOffset=alt.XOffset("region:N", sort=legend_order),
        y=alt.Y(
            "fraction:Q",
            title="Fraction of reads",
            scale=alt.Scale(type="log"),
            stack=None,
            axis=alt.Axis(format="~g"),
        ),
        color=alt.Color(
            "region:N",
            scale=alt.Scale(domain=legend_order, range=["#8D9279", "#05AA8F"]),
            legend=alt.Legend(title=None, orient="top"),
        ),
    )
    .properties(width=350, height=250)
)

chart.save(snakemake.output[0], embed_options={"actions": False}, inline=False)
