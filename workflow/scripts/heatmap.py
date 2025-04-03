import altair as alt
import pandas as pd
import numpy as np
import os
import pysam
import cyvcf2

def filter_candidates(
    bcf_file, bam_file, mapq_threshold=60, min_fraction=0.5
):
    """
    Filter variants from a BCF file based on mapping quality (MAPQ) of aligned reads.
    """
    positions_low_mapq = []
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Open the BCF file
    vcf_reader = cyvcf2.VCF(bcf_file)

    for record in vcf_reader:
        chrom = record.CHROM
        pos = record.POS

        # Fetch reads at the given position
        reads = bam.fetch(chrom, pos - 1, pos)

        mapq_values = [read.mapping_quality for read in reads if not read.is_unmapped]
        if not mapq_values:
            continue  # No reads found

        # Calculate the fraction of reads with MAPQ < threshold
        low_mapq_count = sum(1 for q in mapq_values if q < mapq_threshold)
        fraction_low_mapq = low_mapq_count / len(mapq_values)

        if fraction_low_mapq >= min_fraction:
            positions_low_mapq.append(pos)

    return positions_low_mapq




def compute_rmse(df):
    squared_errors = (df["tool_methylation"] - df["true_methylation"]) ** 2
    mean_squared_error = squared_errors.mean()
    return np.sqrt(mean_squared_error)


def plot_meth_vals_heatmap(df, output, tool_name):
    rmse = compute_rmse(df)
    print(df)
    heatmap = (
        alt.Chart(df)
        .transform_bin("x_bin", field="tool_methylation", bin=alt.Bin(maxbins=100))
        .transform_bin("y_bin", field="true_methylation", bin=alt.Bin(maxbins=100))
        .transform_aggregate(count="count()", groupby=["x_bin", "y_bin"])
        .mark_rect()
        .encode(
            x=alt.X(
                "x_bin:O", 
                title="Tool Methylation",
                axis=alt.Axis(labelExpr="datum.value % 10 == 0 ? datum.value : ''")
            ),
            y=alt.Y(
                "y_bin:O", 
                title="True Methylation", 
                sort="descending",
                axis=alt.Axis(labelExpr="datum.value % 10 == 0 ? datum.value : ''")
            ),
            color=alt.Color("count:Q", scale=alt.Scale(scheme="viridis", type="log"), title="Density"),
            tooltip=["count:Q"]
        )
    )
    line = (
        alt.Chart(pd.DataFrame({"x": [0, 100], "y": [0, 100]}))
        .mark_line(color="red")
        .encode(
            x=alt.X("x:Q", axis=None),  # Entfernt Achsenbeschriftung für X
            y=alt.Y("y:Q", axis=None)   # Entfernt Achsenbeschriftung für Y
        )
    )

    chart = (
        (heatmap + line) 
        .properties(
            width=400,
            height=400,
            title=alt.Title(
                f"{tool_name} vs. TrueMeth",
                subtitle=f"Rmse: {rmse}, datapoints: {len(df)}",
            ),
        )
        .interactive()
    )

    chart.save(output, scale_factor=2.0)


alt.data_transformers.enable("vegafusion")
pd.set_option("display.max_columns", None)



tool_file = snakemake.input[0]
base_name = os.path.splitext(os.path.basename(tool_file))[0]
file_name = "Varlociraptor" if base_name == "calls" else base_name



df = pd.read_parquet(tool_file, engine="pyarrow")
# print("Testcase:", base_name, tool_file, df)
# print(df.to_string())

if file_name == "varlo":
    df = df[df["prob_present"] >= float(snakemake.params["prob_pres_threshhold"])]
cov_bin = int(snakemake.params["cov_bin"])
if cov_bin != -1:
    df = df[df["cov_bin"] == cov_bin]

print(snakemake.params["filter_candidates"])

# if snakemake.params["filter_candidates"] == 1:
#     filtered_candidates = filter_candidates(
#         snakemake.input["candidates"], snakemake.input["alignment"]
#     )
#     with open("/projects/koesterlab/benchmark-methylation/varlociraptor-methylation-evaluation/mapq_output.txt", "w") as f:
#         f.write(str(filtered_candidates))
#     print(df, filtered_candidates)

plot_meth_vals_heatmap(df, snakemake.output["plot"], file_name)


with open(snakemake.input["filtered_candidates"], 'r') as file:
    # Die Textdatei einlesen und die Klammern und Leerzeichen entfernen
    position_list = eval(file.read())

df = df[df['position'].isin(position_list)]
plot_meth_vals_heatmap(df, snakemake.output["plot_filtered"], file_name)