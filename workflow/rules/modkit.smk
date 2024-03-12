rule install_modkit:
    output:
        "modkit_installed.flag",
    conda:
        "../envs/modkit.yaml"
    log:
        "../logs/install_modkit.log",
    shell:
        """
        export PATH=$PATH:~/.cargo/bin
        export PATH=$PATH:/homes/aprinz/.cargo/bin
        cargo install --git https://github.com/nanoporetech/modkit.git --tag v0.2.4
        touch {output}
        """


# Needs a fasta with >chr1 instead of >1
rule modkit:
    input:
        installed="modkit_installed.flag",
        alignment="resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup.bam",
        alignment_index="resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup.bam.bai",
        chromosome=expand(
            "resources/chr_chromosome_{chromosome}.fasta",
            chromosome=chromosome_conf["chromosome"],
        ),
    output:
        "results/Nanopore/{protocol}/alignments_CpG.combined.bed",
    conda:
        "../envs/modkit.yaml"
    log:
        "logs/modkit_{protocol}.log",
    shell:
        """
        export PATH=$PATH:~/.cargo/bin
        export PATH=$PATH:/homes/aprinz/.cargo/bin
        modkit pileup {input.alignment} {output} --cpg --ref {input.chromosome} --force-allow-implicit --combine-strands --log-filepath log_modkit.txt
        echo modkit pileup {input.alignment} {output} --cpg --ref {input.chromosome} --force-allow-implicit --combine-strands --log-filepath log_modkit.txt
        """
