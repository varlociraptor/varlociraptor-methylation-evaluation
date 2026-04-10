# https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/

rule bismark_copy_genome:
    input:
        "resources/{chrom}.fasta",
    output:
        "resources/ref_tools/bismark/genome/{platform}/{chrom}.fasta",
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark/bismark_copy_chromosome/{chrom}_{platform}.log",
    shell:
        """
        mkdir -p $(dirname {output}) 2> {log}
        cp {input} {output} 2> {log}
        """




rule bismark_genome_preparation_fa:
    input:
        genome=lambda wildcards: (
            expand(
                "resources/{chrom}.fasta",
                chrom=config["seq_platforms"].get(wildcards.platform),
            )
            if wildcards.sample.startswith("simulated_data")
            else ["resources/genome.fasta"]
        ),
    output:
        bismark_genome_dir=directory("resources/ref_tools/bismark/genome/{platform}/bismark/"),

    log:
        "logs/bismark_genome_preparation/{platform}.log",
    params:
        extra="",  # optional params string
    threads: 4  # bismark_genome_preparation requires least 2 threads and at least --cores 2 from workflow run
    wrapper:
        "v9.4.1/bio/bismark/bismark_genome_preparation"


rule bismark_align:
    input:
        fq_1="resources/{platform}/{sample}/{SRA}/{SRA}_1_trimmed.fastq",
        fq_2="resources/{platform}/{sample}/{SRA}/{SRA}_2_trimmed.fastq",

        bismark_indexes_dir="resources/ref_tools/bismark/genome/{platform}/",
        # ct="resources/ref_tools/bismark/genome/{platform}/Bisulfite_Genome/CT_conversion",
    output:
        bam="resources/ref_tools/bismark/{platform}/bams/{sample}_pe_{SRA}_unsorted.bam",
        report="resources/ref_tools/bismark/{platform}/bams/{sample}_{SRA}_PE_report.txt",
        fq_unmapped_1="results/ref_tools/bismark/{platform}/{sample}/{SRA}_unmapped_reads_1.fq.gz",  # optional: implicitly activates --unmapped
        fq_unmapped_2="results/ref_tools/bismark/{platform}/{sample}/{SRA}_unmapped_reads_2.fq.gz",  # optional: implicitly activates --unmapped
        fq_ambiguous_1="results/ref_tools/bismark/{platform}/{sample}/{SRA}_ambiguous_reads_1.fq.gz",  # optional: implicitly activates --ambiguous
        fq_ambiguous_2="results/ref_tools/bismark/{platform}/{sample}/{SRA}_ambiguous_reads_2.fq.gz"
    log:
        "logs/bismark/bismark_align/{sample}_{SRA}_{platform}.log",
    benchmark:
        repeat("benchmarks/{platform}/bismark/bismark_align_{SRA}/{sample}.bwa.benchmark.txt", 3)
    params:
        extra="",
    threads: 8
    resources:
        mem_mb=16000,
    wrapper:
        "v9.3.0/bio/bismark/bismark"


# merge bam files from different lanes
rule samtools_merge:
    input:
        get_sample_sra_bismark,
    output:
        "resources/ref_tools/bismark/{platform}/bams/{sample}_pe.bam",
    log:
        "logs/bismark/samtools_merge/{sample}_{platform}.log",
    params:
        extra="-n -f",
    benchmark:
        repeat("benchmarks/{platform}/bismark/samtools_merge/{sample}.bwa.benchmark.txt", 3)
    threads: 8
    wrapper:
        "v5.9.0/bio/samtools/merge"


rule samtools_sort:
    input:
        "resources/ref_tools/bismark/{platform}/bams/{sample}_pe.bam",
    output:
        "resources/ref_tools/bismark/{platform}/bams/{sample}_pe_sorted.bam",
    log:
        "logs/bismark/samtools_sort/{sample}_{platform}.log",
    params:
        extra="-m 4G -n",
    threads: 8
    resources:
        mem_mb=16000,
    benchmark:
        repeat("benchmarks/{platform}/bismark/samtools_sort/{sample}.bwa.benchmark.txt", 3)
    wrapper:
        "v5.9.0/bio/samtools/sort"


rule deduplicate_bismark:
    input:
        "resources/ref_tools/bismark/{platform}/bams/{sample}_pe_sorted.bam",
    output:
        bam="resources/ref_tools/bismark/{platform}/dedup/{sample}.deduplicated.bam",
        report="resources/ref_tools/bismark/{platform}/dedup/{sample}.deduplication_report.txt",
    log:
        "logs/bismark/deduplicate_bismark/{sample}_{platform}.log",
    params:
        extra="",  # optional params string
    benchmark:
        repeat("benchmarks/{platform}/bismark/deduplicate_bismark/{sample}.bwa.benchmark.txt", 3)
    resources:
        mem_mb=16000,
    wrapper:
        "v9.3.0/bio/bismark/deduplicate_bismark"


# rule bismark_methylation_extractor:
#     input:
#         bam="resources/ref_tools/bismark/{platform}/dedup/{sample}.deduplicated.bam",
#     output:
#         cov_zero_based="resources/ref_tools/bismark/{platform}/meth/{sample}.deduplicated.bedGraph.gz.bismark.zero.cov",
#         mbias_r1="resources/ref_tools/bismark/{platform}/qc/meth/{sample}.deduplicated.M-bias_R1.png",
#         # Only for PE BAMS:
#         mbias_r2="resources/ref_tools/bismark/{platform}/qc/meth/{sample}.deduplicated.M-bias_R2.png",
#         mbias_report="resources/ref_tools/bismark/{platform}/report/meth/{sample}.deduplicated.M-bias.txt",
#         splitting_report="resources/ref_tools/bismark/{platform}/report/meth/{sample}.deduplicated_splitting_report.txt",
#         # 1-based start, 1-based end ('inclusive') methylation info: % and counts
#         methylome_CpG_cov="resources/ref_tools/bismark/{platform}/meth/cov/{sample}.deduplicated.bismark.cov.gz",
#         # BedGraph with methylation percentage: 0-based start, end exclusive
#         methylome_CpG_mlevel_bedGraph="resources/ref_tools/bismark/{platform}/meth/bedgraph/{sample}.deduplicated.bedGraph.gz",
#         # Primary output files: methylation status at each read cytosine position: (extremely large)
#         read_base_meth_state_cpg="resources/ref_tools/bismark/{platform}/meth/CpG_context_{sample}.deduplicated.txt.gz",
#         # * You could merge CHG, CHH using: --merge_non_CpG
#         read_base_meth_state_chg="resources/ref_tools/bismark/{platform}/meth/CHG_context_{sample}.deduplicated.txt.gz",
#         read_base_meth_state_chh="resources/ref_tools/bismark/{platform}/meth/CHH_context_{sample}.deduplicated.txt.gz",
#         # cytosine_report="resources/ref_tools/bismark/{platform}/report/meth/{sample}.deduplicated.cytosine_report.txt",
#     log:
#         "logs/bismark/bismark_methylation_extractor/{sample}_{platform}.log",
#     params:
#         output_dir="resources/ref_tools/bismark/{platform}/meth",  # optional output dir
#         extra="--gzip --comprehensive --bedGraph --zero_based",  # optional params string
#     benchmark:
#         repeat("benchmarks/{platform}/bismark/bismark_methylation_extractor/{sample}_{platform}.bwa.benchmark.txt", 3)
#     resources:
#         mem_mb=16000,
#     wrapper:
#         "v5.9.0/bio/bismark/bismark_methylation_extractor"


rule bismark_extract:
    input:
        bam="resources/ref_tools/bismark/{platform}/dedup/{sample}.deduplicated.bam",
    output:
        cov_zero_based="resources/ref_tools/bismark/{platform}/meth/{sample}.deduplicated.bedGraph.gz.bismark.zero.cov",
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark_extract/{sample}_{platform}.log",
    benchmark:
        repeat("benchmarks/{platform}/bismark/bismark_methylation_extractor/{sample}.bwa.benchmark.txt", 3)
    resources:
        mem_mb=16000,
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output}) 2> {log}
        bismark_methylation_extractor {input} -o $(dirname {output}) --parallel {threads} --comprehensive --gzip --comprehensive --bedGraph --zero_based 2> {log}
        """


# We need this rule since the --comprehensive option in bismark_methylation_extractor
# does not create the desired bedGraph file with merged positions for forward and reverse read. We merge them manually by comparing to our candidates.
rule bismark_merge_positions:
    input:
        bedgraph="resources/ref_tools/bismark/{platform}/meth/{sample}.deduplicated.bedGraph.gz.bismark.zero.cov",
        candidates=lambda wildcards: expand(
            "resources/{chrom}/candidates.bcf",
            chrom=config["seq_platforms"].get(wildcards.platform),
        ),
        candidates_index=lambda wildcards: expand(
            "resources/{chrom}/candidates.bcf.csi",
            chrom=config["seq_platforms"].get(wildcards.platform),
        ),
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/bismark.bed",
    log:
        "logs/bismark/bismark_merge_positions/{platform}_{sample}_{platform}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/merge_forward_reverse_positions.py"
