rule methylDackel:
    input:
        #TODO: Wenn Probleme versuche genome statt chromosome, aber dann klappts nicht mehr mit simulated
        # chromosome="resources/genome.fasta",
        chromosome=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["simulated_chrom"],
        ),
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        "results/Illumina_pe/{protocol}/result_files/alignments_CpG.bedGraph",
    conda:
        "../envs/methylDackel.yaml"
    log:
        "logs/methylDackel_{protocol}.log",
    params:
        pipeline_path=config["pipeline_path"],
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ".bedGraph"
        ),
    shell:
        """
        OUTDIR=$(dirname {output})/alignments
        MethylDackel extract {input.chromosome} {input.alignment} -o $OUTDIR --mergeContext
        """


rule rename_dackel_output:
    input:
        "results/Illumina_pe/{protocol}/result_files/alignments_CpG.bedGraph",
    output:
        "results/Illumina_pe/{protocol}/result_files/methylDackel.bed",
    shell:
        "mv {input} {output}"
