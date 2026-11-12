import sys

import altair as alt
import numpy as np
import polars as pl

sys.stderr = open(snakemake.log[0], "w")


def phred_to_prob(score: str) -> float:
    return 10 ** (-float(score) / 10)


def parse_afd(afd: str):
    """Turn the AFD field ("θ=phred,...") into a θ grid and the CDF of the posterior density on ]0,1]."""
    # Duplicate θ values collapse in the dict; sorting orders the grid
    grid = {float(theta): float(prob) for theta, prob in (entry.split("=") for entry in afd.split(","))}
    # θ = 0 is the point mass of the ABSENT event, not a density value; it enters via PROB_ABSENT
    sorted_grid = sorted((theta, phred) for theta, phred in grid.items() if theta > 0)
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


def compute_pit(df: pl.DataFrame, rng) -> pl.DataFrame:
    """Compute u = P(θ ≤ θ* | data) for every site, mixing the point mass at θ = 0 (ABSENT) with the AFD on ]0,1]."""
    u_values = []
    for afd, p_abs, p_low, p_high, theta_true in df.select(
        "afd", "prob_absent", "prob_low", "prob_high", "theta_true"
    ).iter_rows():
        parsed = parse_afd(afd)
        total = p_abs + p_low + p_high
        # We attach the pit values to the original dataframe, so we need to keep the same number of rows.
        if total <= 0:
            u_values.append(None)
            continue
        # Posterior mass of the point θ = 0, conditional on the site not being an artifact
        w0 = p_abs / total
        # An AFD without area on ]0,1] is only usable if all mass lies at θ = 0
        if parsed is None and w0 < 1:
            u_values.append(None)
            continue
        if theta_true <= 0:
            # We need to sample from the uniform distribution on [0, w0] to get a uniform distribution of u values for sites with θ* = 0. Else we would get a spike at u = w0.
            u = rng.uniform(0, w0)
        elif parsed is None:
            # If we have no AFD, but the true methylation level is > 0, we assume that the posterior mass is concentrated at θ = 0.
            u = 1.0
        else:
            theta, cdf = parsed
            u = w0 + (1 - w0) * np.interp(theta_true, theta, cdf)
        u_values.append(float(u))
    return df.with_columns(
        pl.Series("u", u_values, dtype=pl.Float64),
    ).filter(pl.col("u").is_not_null())


truth_df = pl.read_csv(snakemake.input.truth, schema_overrides={"chrom": pl.Utf8})
calls_df = parse_calls(snakemake.input.calls)

num_bins = snakemake.params.num_bins
rng = np.random.default_rng(snakemake.params.seed)

df = calls_df.join(truth_df, on=["chrom", "pos"], how="inner").with_columns(
    (pl.col("methylation_level") / 100).alias("theta_true"),
)
df = compute_pit(df, rng)

df.select(
    "chrom",
    "pos",
    "total_coverage",
    "theta_true",
    "u",
).write_csv(snakemake.output.table)

# Histogram of u with absolute counts per bin
edges = np.linspace(0, 1, num_bins + 1)
counts, _ = np.histogram(df["u"].to_numpy(), bins=edges)
n = int(counts.sum())
hist_df = pl.DataFrame(
    {
        "bin_start": edges[:-1],
        "bin_end": edges[1:],
        "count": counts,
    }
)

ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]

bars = (
    alt.Chart()
    .mark_bar(color="#05AA8F", stroke="white", strokeWidth=0.5)
    .encode(
        x=alt.X(
            "bin_start:Q",
            scale=alt.Scale(domain=[0, 1]),
            axis=alt.Axis(values=ticks),
            title="P(θ ≤ θ* | data)",
        ),
        x2="bin_end:Q",
        y=alt.Y("count:Q", scale=alt.Scale(zero=True), title="number of sites"),
        y2=alt.datum(0),
    )
)
# Count per bin expected under perfect calibration (uniform u)
uniform = (
    alt.Chart()
    .mark_rule(strokeDash=[4, 4], color="#888888", strokeWidth=1)
    .encode(y=alt.datum(n / num_bins))
)

chart = alt.layer(bars, uniform, data=hist_df).properties(
    width=300, height=200, title=f"n = {n}"
)
chart.save(snakemake.output.plot)
