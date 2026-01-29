import pandas as pd

# Show 100 rows for pd
pd.set_option("display.max_rows", 300)

df_untreated = pd.read_parquet(snakemake.input[0])
df_emseq = pd.read_parquet(snakemake.input[1])
print(df_emseq[df_emseq["position"] == 44940672])
df_both = pd.read_parquet(snakemake.input[2])

print(df_both.head())

df_merged = df_untreated.merge(
    df_emseq,
    on=["chromosome", "position", "ref", "alt"],
    suffixes=("_untreated", "_emseq"),
)

df_merged = df_merged.merge(
    df_both, on=["chromosome", "position", "ref", "alt"], suffixes=("", "_both")
)

print(df_merged.columns)
df_merged["diff_probs"] = (
    df_merged["prob_present_emseq"] - df_merged["prob_present_untreated"]
)


# Show differences where diff_probs != 0
df_diff = df_merged[df_merged["diff_probs"] >= 1]

# Keep only prob_present cols
df_diff = df_diff[
    [
        "chromosome",
        "position",
        "prob_present_untreated",
        "prob_present_emseq",
        "prob_present",
        "diff_probs",
    ]
]
print(df_diff.head(300))
