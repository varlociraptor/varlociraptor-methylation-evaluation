# Needs a fasta with >chr1 instead of >1
rule modkit_compute_methylation:
    input:
        alignment="resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup.bam",
        alignment_index="resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup.bam.bai",
        chromosome=expand(
            "resources/chr_chromosome_{chromosome}.fasta",
            chromosome=chromosome_by_seq_platform.get("Nanopore"),
        ),
    output:
        "results/Nanopore/{protocol}/result_files/alignments_CpG.combined.bed",
    conda:
        "../envs/modkit.yaml"
    log:
        "logs/modkit/{protocol}/compute_methylation.log",
    resources:
        mem_mb=128000
    shell:
        """
        export PATH=$PATH:~/.cargo/bin 2> {log}
        export PATH=$PATH:/homes/aprinz/.cargo/bin 2> {log}
        modkit pileup {input.alignment} {output} --cpg --ref {input.chromosome} --force-allow-implicit --combine-strands --log-filepath {log} 2> {log}
        """


rule modkit_rename_output:
    input:
        "results/Nanopore/{protocol}/result_files/alignments_CpG.combined.bed",
    output:
        "results/Nanopore/{protocol}/result_files/modkit.bed",
    log:
        "logs/modkit/{protocol}/rename_output.log",
    shell:
        "mv {input} {output} 2> {log}"
