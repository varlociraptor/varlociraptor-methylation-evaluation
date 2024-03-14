rule bsmap:
    input:
        genome="resources/genome.fasta",
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup.bam.bai",
    output:
        "results/ref_tools/bsMap/{protocol}/out.sam",
    conda:
        "../envs/bsmap.yaml"
    log:
        "logs/bsmapz_{protocol}.log",
    threads: 8
    shell:
        """
        bsmap -a {input.alignment} -b {input.alignment} -d {input.genome} -o out.sam -p {threads} -w 100  -v 0.07 -m 50 -x 300
        cp out.sam {output}
        """


rule extract_methylation:
    input:
        genome="resources/genome.fasta",
        bsmap_sam="results/ref_tools/bsMap/{protocol}/out.sam",
    output:
        "results/ref_tools/bsMap/{protocol}/methylation_ratios.txt",
    conda:
        "../envs/bsmap.yaml"
    log:
        "logs/bsmapz_{protocol}.log",
    params:
        chromosome=chromosome_conf["chromosome"],
    shell:
        """
        python ~/Downloads/methratio.py --chr={params.chromosome} --ref={input.genome} --out={output} out.sam -g -x CG
        """


# bsmap -a resources/Illumina_pe/test/alignment_focused_downsampled_dedup.bam -b resources/Illumina_pe/test/alignment_focused_downsampled_dedup.bam -o out.sam -d resources/genome.fasta  -p 1 -w 100  -v 0.07
#  -m 50 -x 300
# python ~/Downloads/methratio.py --chr=J02459.1 --ref=resources/example_paired/genome.fasta --out={output} {input.bsmap_sam} -g -x CG
