# We need the newer methylation extractor from bsmapz, since the original one is outdated
rule meth_extractor:
    output:
        "resources/ref_tools/bsMap/methratio.py",
    log:
        "../logs/meth_extractor.log",
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p mkdir -p $(dirname {output})
        cd $(dirname {output})
        wget https://raw.githubusercontent.com/zyndagj/BSMAPz/master/methratio.py
        """


# We have to rename the alignment file to make it shorter to not get a buffer error
# TODO: does not work with symlinks on HPC: https://github.com/zyndagj/BSMAPz/issues/19 Filenames are too long.
# Maximum allowed number of characters in alignment name: 68
rule bsmap:
    input:
        genome=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        # alignment_renamed=temp("bsmap.bam"),
        out="results/Illumina_pe/{protocol}/result_files/out.sam",
    conda:
        "../envs/bsmap.yaml"
    threads: 8
    shell:
        """
        cp {input.alignment} temp.bam
        bsmap -a temp.bam -b temp.bam -d {input.genome} -o out.sam -p {threads} -w 100  -v 0.07 -m 50 -x 300
        mv out.sam $(dirname {output.out})
        """
        # cp {input.alignment} {output.alignment_renamed}
        # bsmap -a {output.alignment_renamed} -b {output.alignment_renamed} -d {input.genome} -o out.sam -p {threads} -w 100  -v 0.07 -m 50 -x 300
        # mv out.sam $(dirname {output.out})


rule extract_methylation:
    input:
        genome=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        genome_index=expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        bsmap_sam="results/Illumina_pe/{protocol}/result_files/out.sam",
        meth_extractor="resources/ref_tools/bsMap/methratio.py",
    output:
        "results/Illumina_pe/{protocol}/result_files/methylation_ratios.bed",
    conda:
        "../envs/bsmap.yaml"
    log:
        "logs/bsmap_{protocol}.log",
    params:
        chromosome=chromosome_by_platform.get("Illumina_pe"),
    shell:
        """
        python {input.meth_extractor} -c={params.chromosome} --ref={input.genome} --out={output} {input.bsmap_sam} -g -x CG
        """


rule rename_bsmap_output:
    input:
        "results/Illumina_pe/{protocol}/result_files/methylation_ratios.bed",
    output:
        "results/Illumina_pe/{protocol}/result_files/bsMap.bed",
    shell:
        "mv {input} {output}"
