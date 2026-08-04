# Needs a fasta with >chr1 instead of >1
rule modkit_compute_methylation:
    input:
        alignment="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        chromosome=lambda wildcards: expand(
            "resources/{chromosome}.fasta",
            chromosome=chromosome_by_seq_platform.get(wildcards.platform),
        ),
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/modkit.bed",
    log:
        "logs/modkit/modkit_compute_methylation/{platform}_{sample}.log",
    benchmark:
        repeat(
            "benchmarks/{platform}/modkit/modkit/{sample}.bwa.benchmark.txt",
            config["benchmark_repeats"],
        )
    conda:
        "../envs/modkit.yaml"
    threads: 8
    resources:
        mem_mb=16000,
    shell:
        "modkit pileup {input.alignment} {output} --cpg --ref {input.chromosome} --modified-bases 5mC --threads {threads} --combine-strands --log-filepath {log} 2> {log}"
