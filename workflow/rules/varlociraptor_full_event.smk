# Supplementary analysis: Varlociraptor calls with an additional FULL event (θ = 1).
# All outputs live under results/single_sample_full_event, so the regular single_sample results are untouched.


rule varlociraptor_call_full_event:
    input:
        varlo=varlociraptor_binary(),
        preprocess_obs="results/preprocessed/{seq_platform}/{sample}/normal_{scatteritem}.bcf",
        scenario=workflow.source_path("../scenarios/scenario_full.yaml"),
    output:
        "results/single_sample_full_event/{seq_platform}/called/{sample}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor_single/varlociraptor_call_full_event/{seq_platform}_{sample}_{scatteritem}.log",
    params:
        varlo=varlociraptor_cmd,
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{params.varlo} call variants generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output} 2> {log}"


# Replicate heatmaps (MAP) of the regular calls next to the calls with the FULL event
rule plot_full_event_heatmaps:
    input:
        regular=lambda wildcards: expand(
            "results/single_sample/{{seq_platform}}/called/{rep}/result_files/varlo.bed",
            rep=sorted(config["data"][wildcards.seq_platform]),
        ),
        full_event=lambda wildcards: expand(
            "results/single_sample_full_event/{{seq_platform}}/called/{rep}/result_files/varlo.bed",
            rep=sorted(config["data"][wildcards.seq_platform]),
        ),
    output:
        report(
            "results/single_sample_full_event/{seq_platform}/plots/full_event_heatmap.{plot_type}",
            category="supplementary",
            subcategory=lambda wildcards: seq_platform_to_name[wildcards.seq_platform],
            labels={
                "file_type": "full_event_heatmap",
                "sample": "dummy",
            },
            caption="../report/full_event_heatmap.rst",
        ),
    log:
        "logs/plot_results/plot_full_event_heatmaps/{seq_platform}_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=8000,
    params:
        bin_size=config["heatmap_bin_size"],
        fdr_levels=config["fdr_alpha"],
        platform=lambda wildcards: seq_platform_to_name[wildcards.seq_platform],
    script:
        "../scripts/plot_full_event_heatmaps.py"
