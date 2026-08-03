
rule bwameth_index:
    input:
        "resources/{genome}.fasta",
    output:
        temp(
            multiext(
                "resources/{genome}.fasta.bwameth",
                    ".c2t",
                    ".c2t.amb",
                    ".c2t.ann",
                    ".c2t.bwt",
                    ".c2t.pac",
                    ".c2t.sa",
            ),
        ),
    cache: True
    log:
        "logs/bwameth/bwameth_index/{genome}.log",
    wrapper:
        "v7.3.0/bio/bwameth/index"


rule align_reads_pe:
    input:
        ref=lambda wildcards: (
            expand(
                "resources/{chrom}.fasta",
                chrom=config["seq_platforms"].get(wildcards.platform),
            )
            if wildcards.sample.startswith("simulated_data")
            else ["resources/genome.fasta"]
        ),
        idx=lambda wildcards: (
            multiext(expand(

                    "resources/{chrom}.fasta.bwameth", chrom=config["seq_platforms"].get(wildcards.platform))[0],
            ".c2t",
            ".c2t.amb",
            ".c2t.ann",
            ".c2t.bwt",
            ".c2t.pac",
            ".c2t.sa",
            )
            if wildcards.sample.startswith("simulated_data")
            else multiext(
                "resources/genome.fasta.bwameth",
                ".c2t",
                ".c2t.amb",
                ".c2t.ann",
                ".c2t.bwt",
                ".c2t.pac",
                ".c2t.sa",
            )
        ),
        fq1="resources/{platform}/{sample}/{SRA}/{SRA}_1_trimmed.fastq.gz",
        fq2="resources/{platform}/{sample}/{SRA}/{SRA}_2_trimmed.fastq.gz",
    output:
        "resources/{platform}/{sample}/{SRA}/alignment.bam",
    log:
        "logs/bwameth/align_reads_pe/{platform}_{sample}_{SRA}.log",
    threads: 16
    resources:
        mem_mb=48000
    wrapper:
        "v9.4.1/bio/bwameth/memx"



rule aligned_reads_sort:
    input:
        "resources/{seq_platform}/{sample}/{SRA}/alignment.bam",
    output:
        temp("resources/{seq_platform}/{sample}/{SRA}/alignment_sorted.bam"),
    log:
        "logs/bwameth/align_reads_sort/{seq_platform}_{sample}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools sort -@ {threads}  {input} -o {output} 2> {log}"


rule aligned_reads_index:
    input:
        "resources/{seq_platform}/{sample}/{SRA}/alignment_sorted.bam",
    output:
        temp("resources/{seq_platform}/{sample}/{SRA}/alignment_sorted.bam.bai"),
    log:
        "logs/bwameth/aligned_reads_index/{seq_platform}_{sample}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


rule aligned_reads_focus_on_chromosome:
    input:
        bam="resources/{seq_platform}/{sample}/{SRA}/alignment_sorted.bam",
        index="resources/{seq_platform}/{sample}/{SRA}/alignment_sorted.bam.bai",
    output:
        bam="resources/{seq_platform}/{sample}/{SRA}/alignment_focused.bam",
    log:
        "logs/bwameth/aligned_reads_focus_on_chromosome/{seq_platform}_{sample}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        chromosome=lambda wildcards: (
            f"chr{chromosome_by_seq_platform[wildcards.seq_platform]}"
            if wildcards.seq_platform == "PacBio"
            or wildcards.seq_platform == "Nanopore"
            else "21"
            if chromosome_by_seq_platform[wildcards.seq_platform] == "genome"
            else chromosome_by_seq_platform[wildcards.seq_platform]
        ),
        # whole_genome=lambda wildcards: wildcards.seq_platform == "genome",
    threads: 4
    shell:
        # if [ {params.whole_genome} == True ]; then
        #     samtools view -h -@ {threads} -b -o {output.bam} {input.bam} 2> {log}
        # else
        """
        samtools view -h -@ {threads} -b -o {output.bam} {input.bam} {params.chromosome} 2> {log}
        """
        # fi


rule aligned_reads_markduplicates:
    input:
        bams="resources/{seq_platform}/{sample}/{SRA}/alignment_focused.bam",
    output:
        bam="resources/{seq_platform}/{sample}/{SRA}/alignment_focused_dedup.bam",
        metrics="resources/{seq_platform}/{sample}/{SRA}/alignment_focused_dedup.metrics.txt",
    log:
        "logs/bwameth/aligned_reads_markduplicates/{seq_platform}_{sample}_{SRA}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.6.0/bio/picard/markduplicates"


rule aligned_reads_merge_sras:
    input:
        get_sample_sra,
    output:
        "resources/{seq_platform, [^/]+}/{sample,[^/]+}/alignment_focused_dedup.bam",
    log:
        "logs/bwameth/aligned_reads_merge_sras/{seq_platform}_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    wildcard_constraints:
        seq_platform="(?!multi_sample).*",
    shell:
        "samtools merge {output} {input} 2> {log}"


rule aligned_reads_downsample:
    input:
        "resources/{seq_platform}/{sample}/alignment_focused_dedup.bam",
    output:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup.bam",
    log:
        "logs/bwameth/aligned_reads_downsample/{seq_platform}_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools view -h -@ {threads} -s 0.99 -b -o {output} {input} 2> {log}"


rule aligned_reads_downsampled_index:
    input:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup.bam",
    output:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup.bam.bai",
    log:
        "logs/bwameth/aligned_reads_downsampled_index/{seq_platform}_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


rule aligned_reads_rename_chromosomes:
    input:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup.bam",
    output:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
    log:
        "logs/bwameth/aligned_reads_rename_chromosomes/{seq_platform}_{sample}.log",
    # wildcard_constraints:
        # sample="(?!simulated_data).*",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/rename_chrom_in_bam.py"


rule aligned_reads_renamed_index:
    input:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
    output:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    log:
        "logs/bwameth/aligned_reads_renamed_index/{seq_platform}_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"



rule scatter_candidates_to_bed:
    input:
        "resources/{platform}/candidates_{scatteritem}.bcf",
    output:
        "resources/{platform}/candidates_{scatteritem}.bed",
    log:
        "logs/varlociraptor/scatter_candidates_to_bed/{platform}_{scatteritem}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/candidates_to_bed.py"

rule scatter_aligned_reads:
    input:
        alignment="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        index="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        candidate=lambda wildcards: f"resources/{chromosome_by_seq_platform.get(wildcards.seq_platform) if chromosome_by_seq_platform.get(wildcards.seq_platform) != 'genome' else '21'}/candidates_{wildcards.scatteritem}.bed",
    output:
        "resources/{seq_platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam",
    log:
        "logs/bwameth/aligned_reads_candidates_region/{seq_platform}_{sample}_{scatteritem}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -b -L {input.candidate} {input.alignment} > {output} 2> {log}
        # If there is no overlap use he last read since varlo does not work with no reads.
        if [ $(samtools view -c {output}) -eq 0 ]; then
            samtools view -H {input.alignment} > {output}.temp.sam
            samtools view {input.alignment} | tail -n 1 >> {output}.temp.sam
            samtools view -bS {output}.temp.sam > {output}
            rm {output}.temp.sam
        fi
        """


rule aligned_reads_candidates_region_index:
    input:
        "resources/{seq_platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam",
    output:
        "resources/{seq_platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam.bai",
    log:
        "logs/bwameth/aligned_reads_candidates_region_index/{seq_platform}_{sample}_{scatteritem}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"
