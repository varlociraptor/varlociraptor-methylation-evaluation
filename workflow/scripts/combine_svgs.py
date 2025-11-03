import altair as alt
import pickle

platform = snakemake.params["platform"]

with open(snakemake.input[0], "rb") as f:
    chart1 = pickle.load(f)

with open(snakemake.input[1], "rb") as f:
    chart2 = pickle.load(f)

combined = (
    alt.hconcat(chart1, chart2)
        .resolve_scale(color="independent")
        .properties(
            title=alt.TitleParams(
                text=f"{platform.capitalize()} data",
                anchor="middle",
                fontSize=18,
                fontWeight="bold"
            )
        )
)

combined.save(snakemake.output[0])
