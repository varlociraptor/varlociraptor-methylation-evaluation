rule methylDackel_compute_meth:
    input:
        genome=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        genome_index=expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        "results/Illumina_pe/{protocol}/result_files/alignments_CpG.bedGraph",
    conda:
        "../envs/methylDackel.yaml"
    log:
        "logs/methylDackel/{protocol}/compute_meth.log",
    params:
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ".bedGraph"
        ),
    shell:
        """
        OUTDIR=$(dirname {output})/alignments 2> {log}
        MethylDackel extract {input.genome} {input.alignment} -o $OUTDIR --mergeContext 2> {log}
        """


rule methylDackel_rename_output:
    input:
        "results/Illumina_pe/{protocol}/result_files/alignments_CpG.bedGraph",
    output:
        "results/Illumina_pe/{protocol}/result_files/methylDackel.bed",
    log:
        "logs/methylDackel/{protocol}/rename_output.log",
    shell:
        "mv {input} {output} 2> {log}"
