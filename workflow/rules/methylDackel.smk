rule methylDackel:
    input:
        genome="resources/genome.fasta",
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        "results/Illumina_pe/{protocol}/alignments_CpG.combined.bed",
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
        MethylDackel extract {input.genome} {input.alignment} -o {params.pipeline_path}/results/Illumina_pe/{wildcards.protocol}/alignments --mergeContext
        mv {params.pipeline_path}{params.prefix} {params.pipeline_path}{output}
        """
