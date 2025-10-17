# TODO: Do I have to hardcode the inputs, since the shell command uses each of them individually?
rule call_methylation_together:
    input:
        varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
        emseq1="results/Illumina_pe/EMSeq_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        emseq2="results/Illumina_pe/EMSeq_HG002_LAB02_{replicate}/normal_{scatteritem}.bcf",
        methylseq="results/Illumina_pe/MethylSeq_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        splat="results/Illumina_pe/SPLAT_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        truemethylbs="results/Illumina_pe/TrueMethylBS_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        truemethylox="results/Illumina_pe/TrueMethylOX_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        trueseq="results/Illumina_pe/TruSeq_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        pacbio="results/PacBio/{replicate}/normal_{scatteritem}.bcf",
        # nanopore="results/Nanopore/{replicate}/normal_{scatteritem}.bcf",
        scenario="resources/scenario_common_pacbio.yaml",
    output:
        "results/multi_sample/{replicate}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/multi_sample/{replicate}/call_methylation_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs emseq1={input.emseq1} emseq2={input.emseq2} methylseq={input.methylseq} splat={input.splat} truemethylbs={input.truemethylbs} truemethylox={input.truemethylox} trueseq={input.trueseq} pacbio={input.pacbio} > {output} 2> {log}"


# Rename parquet to make it eligible for plotting scripts
rule common_meth_calling_df:
    input:
        tools=[],
        varlo="results/multi_sample/{fdr}/{protocol}/result_files/varlo.parquet",
    output:
        protocol_df="results/multi_sample/{fdr}/{protocol}/result_files/protocol_df_{plot_type}.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/multi_sample/{fdr}/{protocol}/common_tool_df_{plot_type}.log",
    params:
        plot_type=config["plot_type"],
    resources:
        mem_mb=16000,
    script:
        "../scripts/common_tool_df.py"


rule compute_correlation_tables_common:
    input:
        samples=lambda wildcards: expand(
            "results/multi_sample/{fdr}/{protocol}/result_files/protocol_df_{{plot_type}}.parquet",
            fdr=wildcards.fdr,
            protocol=config["data"]["multi_sample"],
        ),
    output:
        table="results/multi_sample/{fdr}/plots/replicates_{plot_type}.hd5",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/multi_sample/{fdr}/correlation_tables_{plot_type}.log",
    params:
        # plot_type=config["plot_type"],
        meth_callers=lambda wildcards: config["ref_tools"]["multi_sample"] + ["varlo"],
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
