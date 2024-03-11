rule methylDackel:
    input:
        genome="resources/genome.fasta",
        alignment="resources/Illumina_{read_type}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_{read_type}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        "results/Illumina_{read_type}/{protocol}/alignments_CpG.combined.bed",
    conda:
        "../envs/methylDackel.yaml"
    log:
        "logs/methylDackel_{read_type}_{protocol}.log",
    params:
        pipeline_path=config["pipeline_path"],
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ".bedGraph"
        ),
    shell:
        """
        MethylDackel extract {input.genome} {input.alignment} -o {params.pipeline_path}/results/Illumina_{wildcards.read_type}/{wildcards.protocol}/alignments --mergeContext
        mv {params.pipeline_path}{params.prefix} {params.pipeline_path}{output}
        """


rule pb_CpG_tools:
    input:
        pb_tools_dir=directory("../pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu"),
        alignment="resources/PacBio/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/PacBio/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        chromosome=expand(
            "resources/chromosome_{chromosome}.fasta",
            chromosome=chromosome_conf["chromosome"],
        ),
    output:
        "results/PacBio/{protocol}/alignments_CpG.combined.bed",
    log:
        "logs/pb_CpG_tools_{protocol}.log",
    params:
        base_dir=config["base_dir"],
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ""
        ),
    shell:
        """
        {params.base_dir}pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
        --bam {input.alignment} \
        --modsites-mode reference \
        --ref {input.chromosome} \
        --output-prefix {params.prefix} \
        --pileup-mode count \
        --threads 8
        """
        # """
        # {params.base_dir}pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
        # --bam {input.alignment} \
        # --output-prefix {params.prefix} \
        # --model {params.base_dir}pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
        # --threads 8
        #         --modsites-mode reference \
        # --ref {input.chromosome} \
        # """


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
        # cargo install --git https://github.com/nanoporetech/modkit.git
