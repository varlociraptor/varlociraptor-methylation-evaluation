# We need the newer methylation extractor from bsmapz, since the original one is outdated
rule meth_extractor:
    output:
        "resources/ref_tools/bsMap/methratio.py",
    log:
        "../logs/meth_extractor.log",
    params:
        pipeline_path=config["pipeline_path"],
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p mkdir -p $(dirname {output})
        cd $(dirname {output})
        wget https://raw.githubusercontent.com/zyndagj/BSMAPz/master/methratio.py
        """


# We have to rename the alignment file to make it shorte to not get a buffer error
rule bsmap:
    input:
        genome="resources/genome.fasta",
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup.bam",
        # alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup.bam.bai",
    output:
        alignment_renamed=temp("resources/Illumina_pe/{protocol}/bsmap.bam"),
        out="results/ref_tools/bsMap/{protocol}/out.sam",
    conda:
        "../envs/bsmap.yaml"
    log:
        "logs/bsmapz_{protocol}.log",
    threads: 8
    shell:
        """
        cp {input.alignment} {output.alignment_renamed}
        bsmap -a {output.alignment_renamed} -b {output.alignment_renamed} -d {input.genome} -o out.sam -p {threads} -w 100  -v 0.07 -m 50 -x 300
        mv out.sam $(dirname {output.out})
        """


rule extract_methylation:
    input:
        genome="resources/genome.fasta",
        bsmap_sam="results/ref_tools/bsMap/{protocol}/out.sam",
        meth_extractor="resources/ref_tools/bsMap/methratio.py",
    output:
        "results/ref_tools/bsMap/{protocol}/methylation_ratios.bed",
    conda:
        "../envs/bsmap.yaml"
    log:
        "logs/bsmap_{protocol}.log",
    params:
        chromosome=chromosome_by_platform["Illumina_pe"],
    shell:
        """
        python {input.meth_extractor} -c={params.chromosome} --ref={input.genome} --out={output} {input.bsmap_sam} -g -x CG
        """


# bsmap -a resources/Illumina_pe/test/alignment_focused_downsampled_dedup.bam -b resources/Illumina_pe/test/alignment_focused_downsampled_dedup.bam -o out.sam -d resources/genome.fasta  -p 1 -w 100  -v 0.07
#  -m 50 -x 300
# python ~/Downloads/methratio.py --chr=J02459.1 --ref=resources/example_paired/genome.fasta --out={output} {input.bsmap_sam} -g -x CG


rule rename_bsmap_output:
    input:
        "results/ref_tools/bsMap/{protocol}/methylation_ratios.bed",
    output:
        "results/ref_tools/bsMap/{protocol}/bsMap.bed",
    shell:
        "mv {input} {output}"
