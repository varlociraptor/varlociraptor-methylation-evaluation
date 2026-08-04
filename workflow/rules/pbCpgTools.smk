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
    benchmark:
        repeat(
            "benchmarks/{platform}/pb-CpG-tools/pb-CpG-tools/{sample}.bwa.benchmark.txt",
            config["benchmark_repeats"],
        )
    conda:
        "../envs/pbcpgtools.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input, output: output[0].replace(
            ".combined.bed.gz", ""
        ),
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
