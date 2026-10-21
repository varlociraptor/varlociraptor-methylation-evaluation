import sys

import altair as alt
import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 10)
alt.data_transformers.enable("vegafusion")


# Copy of the Varlociraptor branch of pandas_df_from_meth_output.py
def read_tool_file(filepath: str, alpha: float) -> pd.DataFrame:
    """    Parse tool-specific methylation output into a standardized DataFrame.
    """
    records = []

    with open(filepath, "r") as f:
        for line in f:
            # Skip header/comment lines
            if line.startswith("#") or line.startswith("track"):
                continue

            parts = line.strip().split("\t")

            # -----------------------------
            # Varlociraptor format
            # -----------------------------
            chrom = parts[0]
            position = int(parts[1])
            alt_allele = parts[4]
            info_field = parts[7]
            format_field = parts[8]

            if alt_allele != "<METH>":
                continue

            # Parse sample AF values
            format_fields = format_field.split(":")
            try:
                af_index = format_fields.index("AF")
            except ValueError:
                continue

            sample_afs = []
            for sample_fmt in parts[9:]:
                values = sample_fmt.split(":")
                try:
                    af = float(values[af_index])
                    sample_afs.append(af * 100)
                except (ValueError, IndexError):
                    pass

            meth_rate = sum(sample_afs) / len(sample_afs) if sample_afs else 0

            # Parse INFO fields
            info_dict = dict(
                item.split("=", 1)
                for item in info_field.strip().split(";")
                if "=" in item
            )

            def phred_to_prob(score):
                return 10 ** (-float(score) / 10)

            def is_missing(val):
                return val is None or val == "."

            # Determine probability of methylation
            if is_missing(info_dict.get("PROB_PRESENT")):
                if is_missing(info_dict.get("PROB_HIGH")) or is_missing(
                    info_dict.get("PROB_LOW")
                ):
                    print(
                        f"Missing probability info at {chrom}:{position}",
                        file=sys.stderr,
                    )
                    continue
                prob_present = phred_to_prob(
                    info_dict["PROB_HIGH"]
                ) + phred_to_prob(info_dict["PROB_LOW"])
            else:
                prob_present = phred_to_prob(info_dict["PROB_PRESENT"])

            # FULL (θ = 1) is a methylated event as well
            if not is_missing(info_dict.get("PROB_FULL")):
                prob_present += phred_to_prob(info_dict["PROB_FULL"])

            prob_absent = phred_to_prob(info_dict.get("PROB_ABSENT", 0))
            prob_artifact = phred_to_prob(info_dict.get("PROB_ARTIFACT", 0))

            if max(prob_present, prob_absent + prob_artifact) < (1 - alpha):
                print(
                    f"Low confidence site skipped: {chrom}:{position}",
                    file=sys.stderr,
                )
                continue
            records.append([chrom, position, meth_rate, parts[9]])

    return pd.DataFrame(
        records, columns=["chromosome", "position", "tool_methylation", "format"]
    )


# Copy of common_tool_df.py
def common_tool_df(tool_dfs: dict) -> pd.DataFrame:
    dfs = []

    for tool_name, df in tool_dfs.items():
        df = df.rename(
            columns={
                "tool_methylation": f"{tool_name}_methylation",
                "format": f"{tool_name}_format",
            }
        )
        # Set chromosome and position as index before appending
        df = df.set_index(["chromosome", "position"])
        dfs.append(df)

    # TODO: Find a better solution.
    # Since we split alignment files for Varlo, BSMPAz,and Methyldackel for mor efficient computation we compute multiple meth rates for positions at the beginning and end of the alignment file.
    # A weighted average would be better but for now we just drop duplicate positions
    all_dup_positions = set()
    for df in dfs:
        all_dup_positions.update(df.index[df.index.duplicated(keep=False)])
    print(f"Found {len(all_dup_positions)} duplicate positions", file=sys.stderr)
    dfs = [df[~df.index.isin(all_dup_positions)] for df in dfs]

    return pd.concat(dfs, axis=1, join="outer").reset_index()


# Copy of prepare_plot_df.py
def bin_methylation(series: pd.Series, bin_size: int) -> pd.Series:
    """Round methylation values to nearest bin_size and cast to int."""

    return (np.round(series / bin_size) * bin_size).astype(int)


def compute_replicate_counts(df, bin_size, meth_callers):
    caller_counts = []
    distances = []
    for caller in meth_callers:
        rep1 = f"{caller}_methylation_rep1"
        rep2 = f"{caller}_methylation_rep2"
        temp = df[[rep1, rep2]].dropna()
        if temp.empty:
            continue

        rep1_vals = temp[rep1].to_numpy()
        rep2_vals = temp[rep2].to_numpy()

        denom = np.maximum(rep1_vals, rep2_vals)
        mape = np.where(
            (rep1_vals == 0) & (rep2_vals == 0),
            0,
            np.abs(rep1_vals - rep2_vals) / denom,
        ).mean()

        binary_discordance = np.mean((rep1_vals > 0) != (rep2_vals > 0))

        mae = np.abs(rep1_vals - rep2_vals).mean()

        distances.append(
            {
                "meth_caller": caller,
                "mape": mape,
                "mae": mae,
                "binary_discordance": binary_discordance,
            }
        )

        temp = temp.assign(
            rep1_bin=bin_methylation(temp[rep1], bin_size),
            rep2_bin=bin_methylation(temp[rep2], bin_size),
        )

        counts = (
            pd.crosstab(temp["rep1_bin"], temp["rep2_bin"])
            .stack()
            .reset_index(name="count")
        )

        counts["meth_caller"] = caller

        caller_counts.append(counts)

    counts_df = pd.concat(caller_counts, ignore_index=True)
    distances_df = pd.DataFrame(distances)

    return counts_df, distances_df


# Copy of plot_heatmaps.py
def plot_heatmap(
    df: pd.DataFrame,
    meth_caller: str,
    bin_size: int,
    distances: dict,
    title: str,
    max_count: int,
) -> alt.Chart:
    """Log-scaled heatmap for replicate methylation counts."""
    ticks = list(np.logspace(0, np.log10(max_count), num=5).round().astype(int))
    mape = distances.loc[distances["meth_caller"] == meth_caller, "mape"].iloc[0]
    mae = distances.loc[distances["meth_caller"] == meth_caller, "mae"].iloc[0]
    binary_discordance = distances.loc[
        distances["meth_caller"] == meth_caller, "binary_discordance"
    ].iloc[0]
    number_points = df["count"].sum() // 1000

    heatmap = (
        alt.Chart(
            df,
            title=alt.TitleParams(
                text=title,
                subtitle=[
                    f"N = {number_points}k | Dᵣ = {mape:.2f} | Dₐ = {mae:.2f} | Dₛ = {binary_discordance:.2f}",
                ],
                subtitleFontSize=9,
            ),
        )
        .mark_rect()
        .encode(
            x=alt.X(
                "rep1_bin:O", sort=list(range(0, 101, bin_size)), title="Replicate 1"
            ),
            y=alt.Y(
                "rep2_bin:O", sort=list(range(100, -1, -bin_size)), title="Replicate 2"
            ),
            # Empty bins are coloured black explicitly (log scale is undefined at 0)
            color=alt.condition(
                "datum.count == 0",
                alt.value("black"),
                alt.Color(
                    "count:Q",
                    scale=alt.Scale(
                        type="log", scheme="viridis", domain=[1, max_count]
                    ),
                    legend=alt.Legend(
                        title="Count",
                        orient="right",
                        values=ticks,
                        format=",",
                        tickCount=len(ticks),
                    ),
                ),
            ),
            tooltip=["rep1_bin:O", "rep2_bin:O", "count:Q"],
        )
        .properties(width=200, height=200)
        .interactive()
    )
    return heatmap


bin_size = snakemake.params["bin_size"]
fdr_levels = snakemake.params["fdr_levels"]
platform = snakemake.params["platform"]
meth_callers = [f"varlo_{alpha}" for alpha in fdr_levels]
call_names = {
    "regular": "Varlociraptor",
    "full_event": "Varlociraptor with FULL event",
}

results = {}
for call_type in call_names:
    # One common tool df per replicate (REP01, REP02), as in common_tool_df.py
    df1, df2 = [
        common_tool_df(
            {f"varlo_{alpha}": read_tool_file(path, alpha) for alpha in fdr_levels}
        )
        for path in snakemake.input[call_type]
    ]
    # As in merge_replicates.py
    merged_df = pd.merge(
        df1,
        df2,
        on=["chromosome", "position"],
        how="inner",
        suffixes=("_rep1", "_rep2"),
    )
    results[call_type] = compute_replicate_counts(merged_df, bin_size, meth_callers)

max_count = max(counts_df["count"].max() for counts_df, _ in results.values())

# Rows: FDR levels, columns: regular and FULL event calls
heatmaps = [
    plot_heatmap(
        counts_df[counts_df["meth_caller"] == m],
        m,
        bin_size,
        distances,
        f"{call_names[call_type]}, α = {alpha}",
        max_count,
    )
    for alpha, m in zip(fdr_levels, meth_callers)
    for call_type, (counts_df, distances) in results.items()
]

alt.concat(*heatmaps, columns=len(call_names)).resolve_scale(
    x="independent", y="independent", color="shared"
).properties(title=platform).save(
    snakemake.output[0], embed_options={"actions": False}, inline=False
)
