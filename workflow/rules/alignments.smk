# rule bwameth_index_genome:
#     input:
#         "resources/{genome}.fasta",
#     output:
#         "resources/{genome}.fasta.bwameth.c2t",
#         # "resources/{genome}.fasta.bwameth.c2t.0123",
#         # "resources/{genome}.fasta.bwameth.c2t.amb",
#         # "resources/{genome}.fasta.bwameth.c2t.ann",
#         # "resources/{genome}.fasta.bwameth.c2t.bwt.2bit.64",
#         # "resources/{genome}.fasta.bwameth.c2t.pac",
#     log:
#         "logs/alignment/index_genome_{genome}.log",
#     conda:
#         "../envs/bwa-meth.yaml"
#     shell:
#         "bwameth.py index {input} 2> {log}"


rule bwameth_index:
    input:
        "resources/{genome}.fasta",
        # "genome.fasta",
    output:
        multiext(
            "resources/{genome}.fasta.bwameth",
            ".c2t",
            ".c2t.amb",
            ".c2t.ann",
            ".c2t.bwt.2bit.64",
            ".c2t.pac",
            ".c2t.0123",
        ),
    cache: True  # save space and time with between workflow caching (see docs)
    threads: 1
    log:
        "logs/bwameth_index/{genome}.log",
    wrapper:
        "v6.0.1/bio/bwameth/index"


rule align_reads_pe:
    input:
        fasta_index=multiext(
            "resources/genome.fasta.bwameth",
            ".c2t",
            ".c2t.amb",
            ".c2t.ann",
            ".c2t.bwt.2bit.64",
            ".c2t.pac",
            ".c2t.0123",
        ),
        # fasta_index="resources/genome.fasta.bwameth.c2t",
        fasta="resources/genome.fasta",
        reads1="resources/Illumina_pe/{sample}/{SRA}/{SRA}_1_trimmed.fastq",
        reads2="resources/Illumina_pe/{sample}/{SRA}/{SRA}_2_trimmed.fastq",
    output:
        "resources/Illumina_pe/{sample}/{SRA}/alignment.bam",
    conda:
        "../envs/bwa-meth.yaml"
    log:
        "logs/alignment/{sample}/align_reads_pe_{SRA}.log",
    wildcard_constraints:
        sample="(?!simulated_data).*",
    threads: 30
    resources:
        mem_mb=512,
    shell:
        """
        touch {input.fasta_index}
        bwameth.py --reference {input.fasta} {input.reads1} {input.reads2} -t {threads}  | samtools view -b - > {output} 2> {log}
        """


rule align_reads_se:
    input:
        fasta_index="resources/genome.fasta.bwameth.c2t",
        fasta="resources/genome.fasta",
        reads1="resources/Illumina_se/{sample}/{SRA}/{SRA}_trimmed.fastq",
    output:
        "resources/Illumina_se/{sample}/{SRA}/alignment.bam",
    conda:
        "../envs/bwa-meth.yaml"
    log:
        "logs/alignment/{sample}/align_reads_se_{SRA}.log",
    threads: 30
    resources:
        mem_mb=512,
    shell:
        "bwameth.py --threads {threads} --reference {input.fasta} | samtools view -S -b - > {output} 2> {log}"


rule aligned_reads_sort:
    input:
        "resources/{seq_platform}/{sample}/{SRA}/alignment.bam",
    output:
        "resources/{seq_platform}/{sample}/{SRA}/alignment_sorted.bam",
    log:
        "logs/alignment/{seq_platform}/{sample}/sort_reads_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads}  {input} -o {output} 2> {log}"


rule aligned_reads_index:
    input:
        "resources/{seq_platform}/{sample}/{SRA}/alignment_sorted.bam",
    output:
        "resources/{seq_platform}/{sample}/{SRA}/alignment_sorted.bam.bai",
    log:
        "logs/alignment/{seq_platform}/{sample}/index_sort_reads_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


rule aligned_reads_focus_on_chromosome:
    input:
        bam="resources/{seq_platform}/{sample}/{SRA}/alignment_sorted.bam",
        index="resources/{seq_platform}/{sample}/{SRA}/alignment_sorted.bam.bai",
    output:
        bam="resources/{seq_platform}/{sample}/{SRA}/alignment_focused.bam",
    log:
        "logs/alignment/{seq_platform}/{sample}/focus_on_chromosome_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        chromosome=lambda wildcards: (
            f"chr{chromosome_by_seq_platform[wildcards.seq_platform]}"
            if wildcards.seq_platform == "PacBio"
            or wildcards.seq_platform == "Nanopore"
            else chromosome_by_seq_platform[wildcards.seq_platform]
        ),
    threads: 1
    shell:
        "samtools view -b -o {output.bam} {input} {params.chromosome} 2> {log}"


rule aligned_reads_filter_on_mapq:
    input:
        "resources/{seq_platform}/{sample}/{SRA}/alignment_focused.bam",
    output:
        "resources/{seq_platform}/{sample}/{SRA}/alignment_focused_filtered.bam",
    log:
        "logs/alignment/{seq_platform}/{sample}/filter_on_mapq_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        min_quality=config["min_mapping_quality"],
    threads: 1
    shell:
        "samtools view -q {params.min_quality} -b -o {output} {input} 2> {log}"


rule aligned_reads_markduplicates:
    input:
        bams="resources/{seq_platform}/{sample}/{SRA}/alignment_focused_filtered.bam",
    output:
        bam="resources/{seq_platform}/{sample}/{SRA}/alignment_focused_dedup.bam",
        metrics="resources/{seq_platform}/{sample}/{SRA}/alignment_focused_dedup.metrics.txt",
    log:
        "logs/alignment/{seq_platform}/{sample}/markduplicates_{SRA}.log",
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
        "logs/alignment/{seq_platform}/{sample}/merge_sras.log",
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
        "logs/alignment/{seq_platform}/{sample}/downsample.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -s 0.99 -b -o {output} {input} 2> {log}"


rule aligned_reads_downsampled_index:
    input:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup.bam",
    output:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup.bam.bai",
    log:
        "logs/alignment/{seq_platform}/{sample}/downsample_index.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


rule aligned_reads_rename_chromosomes:
    input:
        "resources/{seq_platform}/{sample, [^(?!simulated_data$)]}/alignment_focused_downsampled_dedup.bam",
    output:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
    log:
        "logs/alignment/{seq_platform}/{sample}/rename_chromosome.log",
    wildcard_constraints:
        sample="(?!simulated_data).*",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/rename_chrom_in_bam.py"


rule aligned_reads_renamed_index:
    input:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
    output:
        "resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    log:
        "logs/alignment/{seq_platform}/{sample}/renamed_index.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


rule aligned_reads_candidates_region:
    input:
        alignment="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        index="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        candidate=lambda wildcards: f"resources/{chromosome_by_seq_platform.get(wildcards.seq_platform)}/candidates_{wildcards.scatteritem}.bcf",
    output:
        "resources/{seq_platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam",
    log:
        "logs/alignment/{seq_platform}/{sample}/candidates_region_{scatteritem}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        chromosome=lambda wildcards: chromosome_by_seq_platform.get(
            wildcards.seq_platform
        ),
    shell:
        """
        set +o pipefail

        mkdir -p $(dirname {output})

        start=$(bcftools query -f '%POS\n' {input.candidate} | head -n1)
        end=$(bcftools query -f '%POS\n' {input.candidate} | tail -n1)
        samtools view -b {input.alignment} "{params.chromosome}:$start-$end" > {output}

        if [ $(samtools view -c {output}) -eq 0 ]; then
            samtools view -H {input.alignment} > temp.sam
            samtools view {input.alignment} | tail -n 1 >> temp.sam
            samtools view -bS temp.sam > {output}
            rm temp.sam
        fi
        """


rule aligned_reads_candidates_region_index:
    input:
        "resources/{seq_platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam",
    output:
        "resources/{seq_platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam.bai",
    log:
        "logs/alignment/{seq_platform}/{sample}/candidates_region_index_{scatteritem}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"
