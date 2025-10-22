# TODO: Do I have to hardcode the inputs, since the shell command uses each of them individually?
rule call_methylation_together_np_pb:
    input:
        varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
        pacbio="results/PacBio/{replicate}/normal_{scatteritem}.bcf",
        nanopore="results/Nanopore/{replicate}/normal_{scatteritem}.bcf",
        scenario="resources/scenario_common_nanopore_pacbio.yaml",
    output:
        "results/multi_sample/np_pb/{replicate}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/multi_sample/{replicate}/call_methylation_{scatteritem}.log",
    benchmark:
        "benchmarks/multi_sample/np_pb/{replicate}_{scatteritem}.bwa.benchmark.txt"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs  pacbio={input.pacbio}  nanopore={input.nanopore} > {output} 2> {log}"


rule call_methylation_together_np_trueOX:
    input:
        varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
        trueOX="results/Illumina_pe/TrueMethylOX_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        nanopore="results/Nanopore/{replicate}/normal_{scatteritem}.bcf",
        scenario="resources/scenario_common_nanopore_trueOX.yaml",
    output:
        "results/common_calls/np_trueOX/{replicate}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/common_calls/np_trueOX/{replicate}/call_methylation_{scatteritem}.log",
    benchmark:
        "benchmarks/common_calls/np_trueOX/{replicate}_{scatteritem}.bwa.benchmark.txt"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs  trueOX={input.trueOX}  nanopore={input.nanopore} > {output} 2> {log}"

rule call_methylation_together_pb_trueOX:
    input:
        varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
        trueOX="results/Illumina_pe/TrueMethylOX_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        pacbio="results/PacBio/{replicate}/normal_{scatteritem}.bcf",
        scenario="resources/scenario_common_pacbio_trueOX.yaml",
    output:
        "results/common_calls/pb_trueOX/{replicate}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/common_calls/pb_trueOX/{replicate}/call_methylation_{scatteritem}.log",
    benchmark:
        "benchmarks/common_calls/pb_trueOX/{replicate}_{scatteritem}.bwa.benchmark.txt"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs  trueOX={input.trueOX}  pacbio={input.pacbio} > {output} 2> {log}"


rule common_meth_calling_df:
    input:
        tools=[],
        varlo="results/multi_sample/{fdr}/{protocol}/result_files/varlo.parquet",
    output:
        protocol_df="results/multi_sample/{fdr}/{protocol}/result_files/protocol_df.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/multi_sample/{fdr}/{protocol}/common_tool_df.log",
    params:
        plot_type=config["plot_type"],
    resources:
        mem_mb=16000,
    script:
        "../scripts/common_tool_df.py"


rule compute_correlation_tables_common:
    input:
        samples=lambda wildcards: expand(
            "results/multi_sample/{fdr}/{protocol}/result_files/protocol_df.parquet",
            fdr=wildcards.fdr,
            protocol=config["data"]["multi_sample"],
        ),
    output:
        table="results/multi_sample/{fdr}/plots/replicates.hd5",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/multi_sample/{fdr}/correlation_tables.log",
    params:
        # plot_type=config["plot_type"],
        meth_callers=lambda wildcards: config["ref_tools"].get("multi_sample", []) + ["varlo"],
        correlation_methods=config["correlation_methods"],
    script:
        "../scripts/compute_correlation_tables.py"


# Computes one common df out of all single method dfs
# rule df_multi_sample:
#     input:
#         rep1="results/multi_sample/{fdr}/REP01/result_files/varlo.parquet",
#         rep2="results/multi_sample/{fdr}/REP02/result_files/varlo.parquet",
#     output:
#         protocol_df="results/multi_sample/{fdr}/REP/result_files/protocol_df_{plot_type}.parquet",
#     conda:
#         "../envs/plot.yaml"
#     log:
#         "logs/plots/multi_sample/{fdr}/common_tool_df_{plot_type}.log",
#     params:
#         plot_type=config["plot_type"],
#     resources:
#         mem_mb=16000,
#     script:
#         "../scripts/common_tool_df.py"
# # Compute common heatmap over all Illumina protocols
# rule heatmap_multi_sample:
#     input:
#         "results/multi_sample/{fdr}/plots/replicates_{plot_type}.hd5",
#     output:
#         report(
#             "results/multi_sample/{fdr}/plots/heatmap_all_protocols.{plot_type}",
#             category="multi_sample",
#             subcategory=lambda wildcards: f"{wildcards.fdr}",
#             labels={
#                 "file": "heatmap",
#                 "protocol": "all protocols",
#             },
#         ),
#     conda:
#         "../envs/plot.yaml"
#     resources:
#         mem_mb=32000,
#     log:
#         "logs/plots/{fdr}/heatmap_replicates_{plot_type}.log",
#     params:
#         meth_callers=["varlo"],
#         # protocol="all",
#         protocol="varlo",
#         bin_size=lambda wildcards: config["heatmap_bin_size"],
#         # correlation_method=config["correlation_method"],
#     script:
#         "../scripts/heatmap_multi_sample.py"
