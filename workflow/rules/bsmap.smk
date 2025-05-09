# We need the newer methylation extractor from bsmapz, since the original one is outdated:
# https://github.com/zyndagj/BSMAPz
rule bsmap_download:
    output:
        "resources/ref_tools/bsMap/methratio.py",
    log:
        "../logs/bsmap/download.log",
    conda:
        "../envs/shell_cmds.yaml"
    shell:
        """
        mkdir -p mkdir -p $(dirname {output}) 2> {log}
        cd $(dirname {output}) 2> {log}
        wget https://raw.githubusercontent.com/zyndagj/BSMAPz/master/methratio.py 2> {log}
        """


# We have to rename the alignment file to make it shorter to not get a buffer error
# TODO: does not work with symlinks on HPC: https://github.com/zyndagj/BSMAPz/issues/19 Filenames are too long.
# Maximum allowed number of characters in alignment name: 68
rule bsmap_compute_meth:
    input:
        genome=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"]["Illumina_pe"],
        ),
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        # alignment_renamed=temp("bsmap.bam"),
        out="results/Illumina_pe/{protocol}/result_files/out.sam",
    conda:
        "../envs/bsmap.yaml"
    log:
        "logs/bsmap/{protocol}/compute_meth.log",
    threads: 8
    shell:
        """
        cp {input.alignment} temp.bam 2> {log}
        bsmap -a temp.bam -b temp.bam -d {input.genome} -o out.sam -p {threads} -w 100  -v 0.07 -m 50 -x 300 2> {log}
        mv out.sam $(dirname {output.out}) 2> {log}
        """


rule bsmap_extract_meth:
    input:
        genome=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"]["Illumina_pe"],
        ),
        genome_index=expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=config["seq_platforms"]["Illumina_pe"],
        ),
        bsmap_sam="results/Illumina_pe/{protocol}/result_files/out.sam",
        meth_extractor="resources/ref_tools/bsMap/methratio.py",
    output:
        "results/Illumina_pe/{protocol}/result_files/methylation_ratios.bed",
    conda:
        "../envs/bsmap.yaml"
    log:
        "logs/bsmap/{protocol}/extract_meth.log",
    params:
        chromosome=chromosome_by_seq_platform.get("Illumina_pe"),
    shell:
        "python {input.meth_extractor} -c={params.chromosome} --ref={input.genome} --out={output} {input.bsmap_sam} -g -x CG 2> {log}"


rule bsmap_rename_output:
    input:
        "results/Illumina_pe/{protocol}/result_files/methylation_ratios.bed",
    output:
        "results/Illumina_pe/{protocol}/result_files/bsMap.bed",
    log:
        "logs/bsmap/{protocol}/rename_output.log",
    shell:
        "mv {input} {output} 2> {log}"
