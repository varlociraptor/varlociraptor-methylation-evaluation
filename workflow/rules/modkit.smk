# Needs a fasta with >chr1 instead of >1
rule modkit_compute_methylation:
    input:
        alignment="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        chromosome=lambda wildcards: expand(
            "resources/chromosome_{chromosome}.fasta",
            chromosome=chromosome_by_seq_platform.get(wildcards.platform),
        ),
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/modkit.bed",
    conda:
        "../envs/modkit.yaml"
    log:
        "logs/modkit/{platform}/{sample}/compute_methylation.log",
    resources:
        mem_mb=128000,
    benchmark:
        "benchmarks/{platform}/modkit/modkit/{sample}.bwa.benchmark.txt"
    shell:
        """
        export PATH=$PATH:~/.cargo/bin 2> {log}
        export PATH=$PATH:/homes/aprinz/.cargo/bin 2> {log}
        modkit pileup {input.alignment} {output} --cpg --ref {input.chromosome} --force-allow-implicit --combine-strands --log-filepath {log} 2> {log}
        """


# rule modkit_rename_output:
#     input:
#         "results/single_sample/{platform}/called/{sample}/result_files/alignments_CpG.combined.bed",
#     output:
#         "results/single_sample/{platform}/called/{sample}/result_files/modkit.bed",
#     log:
#         "logs/modkit/{platform}/{sample}/rename_output.log",
#     shell:
#         "mv {input} {output} 2> {log}"
