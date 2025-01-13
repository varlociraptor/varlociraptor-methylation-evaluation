# Needs a fasta with >chr1 instead of >1
rule modkit:
    input:
        alignment="resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup.bam",
        alignment_index="resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup.bam.bai",
        chromosome=expand(
            "resources/chr_chromosome_{chromosome}.fasta",
            chromosome=chromosome_by_platform.get("Nanopore"),
        ),
    output:
        "results/Nanopore/{protocol}/result_files/alignments_CpG.combined.bed",
    conda:
        "../envs/modkit.yaml"
    log:
        "logs/modkit_{protocol}.log",
    shell:
        """
        export PATH=$PATH:~/.cargo/bin
        export PATH=$PATH:/homes/aprinz/.cargo/bin
        modkit pileup {input.alignment} {output} --cpg --ref {input.chromosome} --force-allow-implicit --combine-strands --log-filepath log_modkit.txt
        """


rule rename_modkit_output:
    input:
        "results/Nanopore/{protocol}/result_files/alignments_CpG.combined.bed",
    output:
        "results/Nanopore/{protocol}/result_files/modkit.bed",
    shell:
        "mv {input} {output}"
