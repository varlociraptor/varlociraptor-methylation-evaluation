import altair as alt
import pickle
import pandas as pd


platform = snakemake.params["platform"]

df1 = pd.read_parquet(snakemake.input[0])
df2 = pd.read_parquet(snakemake.input[1])
df_summary = pd.concat([df1, df2], ignore_index=True)
if platform != "Illumina_pe":
    heatmap = (
        alt.Chart(
            df,
            title=alt.Title(
                title,
                subtitle=f"N = {df['count'].sum()} | D = {mapes.get(meth_caller,0)}",
            ),
        )
        .mark_rect()
        .encode(
            x=alt.X(
                "rep1_bin:O",
                sort=list(range(0, 101, bin_size)),
                title="replicate 1",
            ),
            y=alt.Y(
                "rep2_bin:O",
                sort=list(range(100, -1, -bin_size)),
                title="replicate 2",
            ),
            color=alt.Color(
                "count:Q",
                scale=alt.Scale(type="log", scheme="viridis", domain=[1, max_count]),
                legend=alt.Legend(
                    title="Count", orient="right", values=ticks, format=","
                ),
            ),
            tooltip=["rep1_bin", "rep2_bin", "count"],
        )
        .properties(width=200, height=200)
    )
else:
    meth_caller_order = [
        "BSMAPz",
        "Bismark",
        "MethylDackel",
        "Varlociraptor α = 1.0",
        "Varlociraptor α = 0.01",
    ]
    colorblind_safe_palette = [
        "#D81B60",
        "#1E88E5",
        "#FFC107",
        "#D35892",
        "#AC3FE6",
    ]

    bars = (
        alt.Chart(df_summary)
        .mark_bar()
        .encode(
            x=alt.X(
                "sample:N",
                axis=alt.Axis(labelAngle=-30),
                title="Sample",
            ),
            xOffset=alt.XOffset("meth_caller:N", sort=meth_caller_order),
            y=alt.Y("distance:Q", title="Discordance"),
            color=alt.Color(
                "meth_caller:N",
                title="Methylation caller",
                scale=alt.Scale(range=colorblind_safe_palette),
                sort=meth_caller_order,
            ),
        )
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
        )
    )

    combined = bars + labels

combined.save(snakemake.output[0])
