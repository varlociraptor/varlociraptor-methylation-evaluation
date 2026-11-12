import re
import sys
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd

# Write errors and print output into the Snakemake log file
sys.stderr = open(snakemake.log[0], "w")
# vegafusion lets Altair handle large tables (data is aggregated before plotting)
alt.data_transformers.enable("vegafusion")

# Columns that together identify one CpG locus
LOCUS = ["sample", "chromosome", "position"]
# FORMAT fields of the Varlociraptor calls (branch fdr-filter-bins). The parquet only stores the values.
VARLO_FORMAT = ["DP", "AF", "SAOBS", "SROBS", "OBS", "OOBS", "SB", "ROB", "RPB", "SCB", "HE", "ALB", "FEB", "AFD"]
BIAS_LABELS = {
    "SB": "Strand Bias",
    "ROB": "Read Orientation Bias",
    "RPB": "Read Position Bias",
    "SCB": "Soft‑clipped Bias",
    "HE": "Haplotype Error",
    "ALB": "Alt Locus Bias",
    "FEB": "Fragment End Bias",
}
BIAS_COLORS = {
    "Strand Bias": "#D81B60",
    "Alt Locus Bias": "#1E88E5",
    "Read Orientation Bias": "#FFC107",
    "Read Position Bias": "#004D40",
    "Soft‑clipped Bias": "#8E44AD",
    "Haplotype Error": "#FF7F0E",
    "Fragment End Bias": "#6D4C41",
}
# Bin edges for the allele frequency (AF) of the other replicate
AF_EDGES = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
# Bar labels, e.g. "Bias, AF ]0.0, 0.2]" (zip pairs each edge with the next one)
CATEGORIES = ["Bias both reps", "Bias, AF = 0"] + [f"Bias, AF ]{lo}, {hi}]" for lo, hi in zip(AF_EDGES, AF_EDGES[1:])]


def artifact_file_info(path):
    """(sample, replicate number, fdr) of .../called/{replicate}/result_files/varlo_artifact_{fdr}.parquet.
    Replicate names without sample prefix belong to "dummy"."""
    path = Path(path)
    # Split the replicate folder (parts[-3]), e.g. "TruSeq_HG002_LAB01_REP01", into sample and number
    m = re.fullmatch(r"(?:(.+)_)?REP0*(\d+)", path.parts[-3])
    # group(1) is None without sample prefix; stem is the file name without ".parquet"
    return m.group(1) or "dummy", int(m.group(2)), path.stem.removeprefix("varlo_artifact_")



def bias_dataframe(fdr):
    """Loci with a confident bias in at least one replicate: passing the ARTIFACT FDR control at level fdr and having
    a bias field (at α = 1.0 the control lets every record through, so only the bias field decides)."""
    # Copy, so the shared table base stays unchanged
    df = base.copy()
    # Loci passing the ARTIFACT FDR control, collected per replicate number
    loci = {1: [], 2: []}
    for path in snakemake.input["artifact"]:
        sample, rep, level = artifact_file_info(path)
        if level == fdr and sample in samples:
            # The files only store chromosome and position, so add the sample column
            loci[rep].append(pd.read_parquet(path, columns=["chromosome", "position"]).assign(sample=sample)[LOCUS])
    for rep, dfs in loci.items():
        # True for each row of base whose locus passed (False if there is no file for this level)
        passed = base_loci.isin(pd.MultiIndex.from_frame(pd.concat(dfs))) if dfs else False
        # Confident = passed the control and has a bias field
        df[f"rep{rep}_confident"] = passed & df[f"bias_rep{rep}"].notna()
    # Keep loci with a confident bias in at least one replicate
    df = df[df["rep1_confident"] | df["rep2_confident"]]
    r1, r2 = df["rep1_confident"], df["rep2_confident"]
    # AF of the replicate without a confident bias
    af_other = np.where(r1, df["AF_rep2"], df["AF_rep1"])
    df = df.assign(
        # Bias type of replicate 1 if it is confident, else of replicate 2
        bias_type_label=df["bias_rep1"].where(r1, df["bias_rep2"]),
        # np.select gives each row the label of its first true condition
        category=np.select(
            [r1 & r2, af_other == 0] + [(af_other > lo) & (af_other <= hi) for lo, hi in zip(AF_EDGES, AF_EDGES[1:])],
            CATEGORIES,
            default="No bias",
        ),
    )
    return df[LOCUS + ["AF_rep1", "AF_rep2", "category", "bias_type_label"]]


def bias_chart(df, fdr):
    if df.empty:
        return alt.Chart(pd.DataFrame({"msg": [f"No bias data for ARTIFACT FDR α = {fdr}"]})).mark_text(size=20).encode(text="msg:N")
    present = [t for t in BIAS_COLORS if t in set(df["bias_type_label"])]
    return (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X("category:N", axis=alt.Axis(labelAngle=-45), title=None, scale=alt.Scale(domain=CATEGORIES)),
            y=alt.Y("count():Q", title="Number of loci"),
            color=alt.Color(
                "bias_type_label:N",
                scale=alt.Scale(domain=present, range=[BIAS_COLORS[t] for t in present]),
                title="Bias Type",
            ),
            tooltip=["category", "count()", "bias_type_label"],
        )
        .properties(width=150, height=200, title=f"ARTIFACT FDR α = {fdr}")
    )


samples = snakemake.params["sample"]
samples = [samples] if isinstance(samples, str) else samples
fmt = ["varlo_1.0_format_rep1", "varlo_1.0_format_rep2"]
table = pd.read_parquet(snakemake.input["table"], columns=LOCUS + fmt)
# Keep the selected samples and the loci called in both replicates
table = table[table["sample"].isin(samples)].dropna(subset=fmt).reset_index(drop=True)
base = table[LOCUS].copy()
for rep in (1, 2):
    # Each FORMAT value is a string like "12:0.5:..."; split it into one column per field
    fields = table[f"varlo_1.0_format_rep{rep}"].str.split(":", expand=True)
    fields.columns = VARLO_FORMAT  # fails if the number of FORMAT fields changed
    # True where a bias field is set ("." means not set)
    biases = fields[list(BIAS_LABELS)].replace(".", np.nan).notna()
    base[f"AF_rep{rep}"] = fields["AF"].astype(float)
    # idxmax: first set bias column (Varlociraptor sets at most one per replicate); where: NaN if none is set
    base[f"bias_rep{rep}"] = biases.idxmax(axis=1).map(BIAS_LABELS).where(biases.any(axis=1))
# (sample, chromosome, position) tuples for fast lookups with isin in bias_dataframe
base_loci = pd.MultiIndex.from_frame(base[LOCUS])
bias_dfs = {fdr: bias_dataframe(fdr) for fdr in snakemake.params["fdrs"]}
platform = snakemake.params["platform"]
alt.hconcat(*[bias_chart(df, fdr) for fdr, df in bias_dfs.items()]).properties(
    title="Illumina" if platform == "Illumina_pe" else platform
).save(snakemake.output[0])
pd.concat([df.assign(fdr=fdr) for fdr, df in bias_dfs.items()], ignore_index=True).to_parquet(snakemake.output[1])
