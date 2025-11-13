import altair as alt
import pickle
import pandas as pd


plot_type = snakemake.params["plot_type"]

print(plot_type)
if plot_type != "bias":
    print("TEST")

    platform = snakemake.params["platform"]

    with open(snakemake.input[0], "rb") as f:
        chart1 = pickle.load(f)

    with open(snakemake.input[1], "rb") as f:
        chart2 = pickle.load(f)
    if platform == "Illumina_pe":
        combined = alt.hconcat(chart1, chart2).resolve_scale(color="independent")
    else:
        combined = (
            alt.hconcat(chart1, chart2)
            .resolve_scale(color="independent")
            .properties(
                title=alt.TitleParams(
                    text=f"{platform.capitalize()} data",
                    anchor="middle",
                    fontSize=18,
                    fontWeight="bold",
                )
            )
        )

    combined.save(snakemake.output[0])
else:
    df1 = pd.read_parquet(snakemake.input[0])
    df2 = pd.read_parquet(snakemake.input[1])
    df_summary = pd.concat([df1, df2], ignore_index=True)
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
