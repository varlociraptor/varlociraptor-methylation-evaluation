# TODO: Do I have to hardcode the inputs, since the shell command uses each of them individually?
rule call_methylation_together_np_pb:
    input:
        pacbio="results/preprocessed/PacBio/{replicate}/normal_{scatteritem}.bcf",
        nanopore="results/preprocessed/Nanopore/{replicate}/normal_{scatteritem}.bcf",
        scenario="resources/scenarios/scenario_common_nanopore_pacbio.yaml",
    output:
        "results/multi_sample/np_pb/called/{replicate}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor_multi/call_methylation_together_np_pb/{replicate}_{scatteritem}.log",
    benchmark:
        "benchmarks/multi_sample/np_bp/np_pb/{replicate}_{scatteritem}.bwa.benchmark.txt"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor call variants generic --scenario {input.scenario} --obs  pacbio={input.pacbio}  nanopore={input.nanopore} > {output} 2> {log}"


rule call_methylation_together_np_trueOX:
    input:
        trueOX="results/preprocessed/Illumina_pe/TrueMethylOX_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        nanopore="results/preprocessed/Nanopore/{replicate}/normal_{scatteritem}.bcf",
        scenario="resources/scenarios/scenario_common_nanopore_trueOX.yaml",
    output:
        "results/multi_sample/np_trueOX/called/{replicate}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor_multi/call_methylation_together_np_trueOX/{replicate}_{scatteritem}.log",
    benchmark:
        "benchmarks/multi_sample/np_trueOX/np_trueOX/{replicate}_{scatteritem}.bwa.benchmark.txt"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor call variants generic --scenario {input.scenario} --obs  trueOX={input.trueOX}  nanopore={input.nanopore} > {output} 2> {log}"


rule call_methylation_together_pb_trueOX:
    input:
        trueOX="results/preprocessed/Illumina_pe/TrueMethylOX_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        pacbio="results/preprocessed/PacBio/{replicate}/normal_{scatteritem}.bcf",
        scenario="resources/scenarios/scenario_common_pacbio_trueOX.yaml",
    output:
        "results/multi_sample/pb_trueOX/called/{replicate}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor_multi/call_methylation_together_pb_trueOX/{replicate}_{scatteritem}.log",
    benchmark:
        "benchmarks/multi_sample/pb_trueOX/pb_trueOX/{replicate}_{scatteritem}.bwa.benchmark.txt"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor call variants generic --scenario {input.scenario} --obs  trueOX={input.trueOX}  pacbio={input.pacbio} > {output} 2> {log}"


rule common_meth_calling_df:
    input:
        tools=[],
        varlo="results/multi_sample/{fdr}/{sample}/result_files/varlo.parquet",
    output:
        sample_df="results/multi_sample/{fdr}/{sample}/result_files/sample_df.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/common_tool_df/{fdr}_{sample}.log",
    params:
        plot_type=config["plot_type"],
    resources:
        mem_mb=16000,
    script:
        "../scripts/common_tool_df.py"


rule compute_correlation_tables_common:
    input:
        samples=lambda wildcards: expand(
            "results/multi_sample/{fdr}/{sample}/result_files/sample_df.parquet",
            fdr=wildcards.fdr,
            sample=config["data"]["multi_sample"],
        ),
    output:
        table="results/multi_sample/{fdr}/plots/replicates.hd5",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/correlation_tables/{fdr}.log",
    params:
        # plot_type=config["plot_type"],
        meth_callers=lambda wildcards: config["ref_tools"].get("multi_sample", [])
        + ["varlo"],
        correlation_methods=config["correlation_methods"],
    script:
        "../scripts/compute_correlation.py"
