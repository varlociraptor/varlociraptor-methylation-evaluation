import altair as alt
import pickle
import pandas as pd


def plot_heatmap(
    df: pd.DataFrame,
    meth_caller: str,
    bin_size: int,
    mapes: dict,
    meth_caller_name: str,
) -> alt.Chart:
    """Log-scaled heatmap for replicate methylation counts."""
    max_count = df["count"].max()
    min_count = 1

    ticks = list(np.logspace(0, np.log10(max_count), num=5).round().astype(int))
    heatmap = (
        alt.Chart(
            df,
            title=alt.Title(
                meth_caller_name,
                subtitle=f"N = {df['count'].sum()} | D = {mapes.get(meth_caller,0):.2f}%",
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
            color=alt.Color(
                "count:Q",
                scale=alt.Scale(type="log", scheme="viridis", domain=[1, max_count]),
                legend=alt.Legend(
                    title="Count",
                    orient="right",
                    values=ticks,
                    format=",",
                    tickCount=len(ticks),
                ),
            ),
            tooltip=["rep1_bin:O", "rep2_bin:O", "count:Q"],
        )
        .properties(width=200, height=200)
        .interactive()
    )
    return heatmap


plot_type = snakemake.params["plot_type"]

if plot_type == "heatmap":
    df1 = pd.read_parquet(snakemake.input[0])
    df2 = pd.read_parquet(snakemake.input[1])
    df = pd.concat([df1, df2], ignore_index=True)
    platform = snakemake.params["platform"]

    combined.save(snakemake.output[0])
elif plot_type == "bias":
    df1 = pd.read_parquet(snakemake.input[0])
    df2 = pd.read_parquet(snakemake.input[1])
    df2 = pd.read_parquet(snakemake.input[2])
    print(df1)
    print(df2)
    print(df2)
    combined = pd.concat([df1, df2, df2], ignore_index=True)

    combined.save(snakemake.output[0])
else:
    df1 = pd.read_parquet(snakemake.input[0])
    df2 = pd.read_parquet(snakemake.input[1])
    df_summary = pd.concat([df1, df2], ignore_index=True)
    meth_caller_order = [
        "Bismark",
        "BSMAPz",
        "MethylDackel",
        "Varlociraptor α = 1.0",
        "Varlociraptor α = 0.01",
    ]
    colorblind_safe_palette = [
        "#D81B60",
        "#1E88E5",
        "#FFC107",
        "#05AA8F",
        "#004D40",
    ]
    bars = (
        alt.Chart(df_summary)
        .mark_bar()
        .encode(
            x=alt.X(
                "sample:N",
                axis=alt.Axis(labelAngle=-30),
                title=None,
            ),
            xOffset=alt.XOffset("meth_caller:N", sort=meth_caller_order),
            y=alt.Y("distance:Q", title="Discordance"),
            color=alt.Color(
                "meth_caller:N",
                title="Methylation caller",
                scale=alt.Scale(range=colorblind_safe_palette),
                sort=meth_caller_order,
            ),
            tooltip=["sample:N", "meth_caller:N", "distance:Q", "number:Q"],
        )
        .interactive()
    )

    labels = (
        alt.Chart(df_summary)
        .mark_text(size=8, dy=-5, color="black")
        .transform_calculate(text_k="datum.number + 'k'")
        .encode(
            text="text_k:N",
            x="sample:N",
            xOffset=alt.XOffset("meth_caller:N", sort=meth_caller_order),
            y="distance:Q",
            tooltip=["sample:N", "meth_caller:N", "distance:Q", "number:Q"],
        )
        .interactive()
    )

    combined = bars + labels

combined.save(snakemake.output[0])
