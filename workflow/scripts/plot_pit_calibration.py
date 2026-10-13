import sys

import altair as alt
import numpy as np
import polars as pl

sys.stderr = open(snakemake.log[0], "w")


def phred_to_prob(score: str) -> float:
    return 10 ** (-float(score) / 10)


def parse_afd(afd: str):
    """Turn the AFD field ("θ=phred,...") into a θ grid and the CDF of the posterior density on it."""
    # Duplicate θ values collapse in the dict; sorting orders the grid
    grid = {float(theta): float(prob) for theta, prob in (entry.split("=") for entry in afd.split(","))}
    sorted_grid = sorted(grid.items())
    if len(sorted_grid) < 2:
        return None

    theta, phred = np.array(sorted_grid).T
    density = 10 ** (-phred / 10)
    # Area between neighbouring grid points (trapezoidal rule), accumulated to the CDF
    areas = np.diff(theta) * (density[1:] + density[:-1]) / 2
    # The CDF at the first grid point is 0
    cdf = np.concatenate([[0.0], np.cumsum(areas)])
    if cdf[-1] <= 0:
        return None
    # Sampled densities do not integrate to 1
    return theta, cdf / cdf[-1]


def parse_calls(filepath: str) -> pl.DataFrame:
    """Extract the AFD and the event probabilities for every call."""
    records = []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[4] != "<METH>":
                continue
            info_dict = dict(
                item.split("=", 1) for item in parts[7].split(";") if "=" in item
            )
            probs = [info_dict.get(f"PROB_{e}") for e in ["ABSENT", "LOW", "HIGH"]]
            if any(p is None or p == "." for p in probs):
                continue
            format_dict = dict(zip(parts[8].split(":"), parts[9].split(":")))
            afd = format_dict.get("AFD")
            if afd is None or afd == ".":
                continue
            p_abs, p_low, p_high = (phred_to_prob(p) for p in probs)
            records.append((parts[0], int(parts[1]), p_abs, p_low, p_high, afd))
    return pl.DataFrame(
        records,
        schema=["chrom", "pos", "prob_absent", "prob_low", "prob_high", "afd"],
    )


def coverage_label(bins, i):
    if i == len(bins) - 1:
        return f"≥ {bins[i]}"
    return f"{bins[i]}–{bins[i + 1] - 1}"


def compute_pit(df: pl.DataFrame, rng) -> pl.DataFrame:
    """Compute u = P(θ ≤ θ* | data) for every site, mixing the point mass at θ = 0 (ABSENT) with the AFD on ]0,1]."""
    u_values = []
    for afd, p_abs, p_low, p_high, theta_true in df.select(
        "afd", "prob_absent", "prob_low", "prob_high", "theta_true"
    ).iter_rows():
        parsed = parse_afd(afd)
        total = p_abs + p_low + p_high
        # We attach the pit values to the original dataframe, so we need to keep the same number of rows.
        if parsed is None or total <= 0:
            u_values.append(None)
            continue
        theta, cdf = parsed
        # Posterior mass of the point θ = 0, conditional on the site not being an artifact
        w0 = p_abs / total
        if theta_true <= 0:
            # We need to sample from the uniform distribution on [0, w0] to get a uniform distribution of u values for sites with θ* = 0. Else we would get a spike at u = w0.
            u = rng.uniform(0, w0)
        else:
            u = w0 + (1 - w0) * np.interp(theta_true, theta, cdf)
        u_values.append(float(u))
    return df.with_columns(
        pl.Series("u", u_values, dtype=pl.Float64),
    ).filter(pl.col("u").is_not_null())


truth_df = pl.read_csv(snakemake.input.truth, schema_overrides={"chrom": pl.Utf8})
calls_df = parse_calls(snakemake.input.calls)

coverage_bins = snakemake.params.coverage_bins
num_bins = snakemake.params.num_bins
rng = np.random.default_rng(snakemake.params.seed)

coverage_labels = [coverage_label(coverage_bins, i) for i in range(len(coverage_bins))]
df = calls_df.join(truth_df, on=["chrom", "pos"], how="inner").with_columns(
    (pl.col("methylation_level") / 100).alias("theta_true"),
    # Assign a coverage bin to each single coverage value
    pl.col("total_coverage")
    .cut(coverage_bins[1:], labels=coverage_labels, left_closed=True)
    .cast(pl.Utf8)
    .alias("coverage"),
)
df = compute_pit(df, rng)

df.select(
    "chrom",
    "pos",
    "coverage",
    "total_coverage",
    "theta_true",
    "u",
).write_csv(snakemake.output.table)

# Histogram of u per coverage bin, normalized to a density (uniform = 1)
strata = [
    (label, df.filter(pl.col("coverage") == label)) for label in coverage_labels
]
edges = np.linspace(0, 1, num_bins + 1)
rows = []
for stratum, stratum_df in strata:
    # For each stratum, compute the histogram of u values and normalize to a density
    counts, _ = np.histogram(stratum_df["u"].to_numpy(), bins=edges)
    n = counts.sum()
    for i, count in enumerate(counts):
        rows.append(
            {
                "coverage": f"{stratum} (n={n})",
                "bin_start": edges[i],
                "bin_end": edges[i + 1],
                "density": count / (n * (edges[i + 1] - edges[i])) if n else 0.0,
            }
        )
hist_df = pl.DataFrame(rows)
facet_order = hist_df["coverage"].unique(maintain_order=True).to_list()

ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]

bars = (
    alt.Chart()
    .mark_bar(color="#05AA8F", stroke="white", strokeWidth=0.5)
    .encode(
        x=alt.X(
            "bin_start:Q",
            scale=alt.Scale(domain=[0, 1]),
            axis=alt.Axis(values=ticks),
            title="u = P(θ ≤ θ* | data)",
        ),
        x2="bin_end:Q",
        y=alt.Y("density:Q", title="density"),
    )
)
uniform = (
    alt.Chart()
    .mark_rule(strokeDash=[4, 4], color="#888888", strokeWidth=1)
    .encode(y=alt.datum(1))
)

chart = (
    alt.layer(bars, uniform, data=hist_df)
    .properties(width=200, height=150)
    .facet(
        facet=alt.Facet("coverage:N", sort=facet_order, title="Coverage"),
        columns=4,
    )
)
chart.save(snakemake.output.plot)
