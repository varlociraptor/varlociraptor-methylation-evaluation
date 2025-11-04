# We need the newer methylation extractor from bsmapz, since the original one is outdated:
# https://github.com/zyndagj/BSMAPz
rule bsmap_download:
    output:
        "resources/ref_tools/bsMap/methratio.py",
    log:
        "logs/bsmap/download.log",
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p logs/bsmap
        cd $(dirname {output})
        wget https://raw.githubusercontent.com/zyndagj/BSMAPz/master/methratio.py 
        """


rule bsmap_compute_meth:
    input:
        genome=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        alignment="resources/Illumina_pe/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        # alignment_renamed=temp("bsmap.bam"),
        temp("results/single_sample/Illumina_pe/called/{sample}/result_files/out.sam"),
    conda:
        "../envs/bsmap.yaml"
    log:
        "logs/bsmap/{sample}/compute_meth.log",
    resources:
        mem_mb=32000,
    benchmark:
        "benchmarks/Illumina_pe/bsmap/bsmap_compute/{sample}.bwa.benchmark.txt"
    threads: 8
    shell:
        """
        bsmap -a {input.alignment} -d {input.genome} -o {output} -p {threads} -w 100  -v 0.07 -m 50 -x 300
        """


rule bsmap_extract:
    input:
        genome=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        genome_index=expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        bsmap_sam="results/single_sample/Illumina_pe/called/{sample}/result_files/out.sam",
        meth_extractor="resources/ref_tools/bsMap/methratio.py",
    output:
        "results/single_sample/Illumina_pe/{sample}/result_files/methylation_ratios.bed",
    conda:
        "../envs/bsmap.yaml"
    log:
        "logs/bsmap/{sample}/extract_meth.log",
    params:
        chromosome=chromosome_by_seq_platform.get("Illumina_pe"),
    benchmark:
        "benchmarks/Illumina_pe/bsmap/bsmap_extract/{sample}.bwa.benchmark.txt"
    shell:
        "python {input.meth_extractor} -c={params.chromosome} --ref={input.genome} --out={output} {input.bsmap_sam} -g -x CG 2> {log}"


rule bsmap_rename_output:
    input:
        "results/single_sample/Illumina_pe/{sample}/result_files/methylation_ratios.bed",
    output:
        "results/single_sample/Illumina_pe/called/{sample}/result_files/bsMap.bed",
    log:
        "logs/bsmap/{sample}/rename_output.log",
    shell:
        """
        mkdir -p $(dirname {output})
        cp {input} {output} 2> {log}
        rm {input} 2> {log}
        """
