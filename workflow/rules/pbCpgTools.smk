# TODO MissingOutputException, even if the output is there. You have to rerun the pipeline afterwards
# rule pb_CpG_download:
#     output:
#         "resources/ref_tools/pb-CpG-tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores",
#         "resources/ref_tools/pb-CpG-tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite",
#     log:
#         "logs/pb_CpG_tools/pb_CpG_download/download.log",
#     conda:
#         "../envs/shell_cmds.yaml"
#     shell:
#         """

#         mkdir -p resources/ref_tools/pb-CpG-tools
#         cd resources/ref_tools/pb-CpG-tools
#         wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.1/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu.tar.gz
#         tar -xzf pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu.tar.gz
#         """


rule pb_CpG_compute_methylation:
    input:
        alignment="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        chromosome=lambda wildcards: expand(
            "resources/{chromosome}.fasta",
            chromosome=chromosome_by_seq_platform.get(wildcards.platform),
        ),
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/alignments_CpG.combined.bed.gz",
    log:
        "logs/pb_CpG_tools/pb_CpG_compute_methylation/{platform}_{sample}.log",
    params:
        prefix=lambda wildcards, input, output: output[0].replace(".combined.bed.gz", "")
    threads: 8
    conda:
        "../envs/pbcpgtools.yaml"
    benchmark:
        repeat("benchmarks/{platform}/pb-CpG-tools/pb-CpG-tools/{sample}.bwa.benchmark.txt", config["benchmark_repeats"])
    shell:
        "aligned_bam_to_cpg_scores --bam {input.alignment} --output-prefix {params.prefix} --threads {threads} 2> {log}"



rule pb_CpG_rename_output:
    input:
        "results/single_sample/{platform}/called/{sample}/result_files/alignments_CpG.combined.bed",
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/pb_CpG_tools.bed",
    log:
        "logs/pb_CpG_tools/pb_CpG_rename_output/{platform}_{sample}.log",
    conda:
        "../envs/general.yaml"
    shell:
        "mv {input} {output} 2> {log}"
