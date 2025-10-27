rule methylDackel_compute_meth:
    input:
        # genome=expand(
        #     "resources/genome.fasta",
        #     chrom=config["seq_platforms"].get("Illumina_pe"),
        # ),
        # genome_index=expand(
        #     "resources/genome.fasta.fai",
        #     chrom=config["seq_platforms"].get("Illumina_pe"),
        # ),
        genome=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        genome_index=expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        alignment="resources/Illumina_pe/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        "results/single_sample/Illumina_pe/called/{sample}/result_files/alignments_CpG.bedGraph",
    conda:
        "../envs/methylDackel.yaml"
    log:
        "logs/methylDackel/{sample}/compute_meth.log",
    benchmark:
        "benchmarks/Illumina_pe/methylDackel/{sample}.bwa.benchmark.txt"
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
        "results/single_sample/Illumina_pe/called/{sample}/result_files/alignments_CpG.bedGraph",
    output:
        "results/single_sample/Illumina_pe/called/{sample}/result_files/methylDackel.bed",
    log:
        "logs/methylDackel/{sample}/rename_output.log",
    shell:
        "mkdir -p $(dirname {output}) && mv {input} {output} 2> {log}"
