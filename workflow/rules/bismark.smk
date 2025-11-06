# https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/


rule bismark_copy_genome:
    input:
        "resources/genome.fasta",
    output:
        "resources/ref_tools/bismark/genome.fasta",
    log:
        "logs/bismark/bismark_copy_genome/copy.log",
    conda:
        "../envs/general.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) 2> {log}
        cp {input} {output} 2> {log}
        """


rule bismark_copy_chromosome:
    input:
        "resources/chromosome_{chrom}.fasta",
    output:
        "resources/ref_tools/bismark/chromosome_{chrom}.fasta",
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark/bismark_copy_chromosome/{chrom}.log",
    shell:
        """
        mkdir -p $(dirname {output}) 2> {log}
        cp {input} {output} 2> {log}
        """


# TODO: Gives missing output exception
rule bismark_prepare_genome:
    input:
        expand(
            "resources/ref_tools/bismark/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
    output:
        directory("resources/ref_tools/bismark/Bisulfite_Genome"),
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark/prepare_genome/prepare.log",
    shell:
        """
        bismark_genome_preparation $(dirname {input}) --verbose  2> {log}
        """


rule bismark_align:
    input:
        fq_1="resources/Illumina_pe/{sample}/{SRA}/{SRA}_1.fastq",
        fq_2="resources/Illumina_pe/{sample}/{SRA}/{SRA}_2.fastq",
        genome=expand(
            "resources/ref_tools/bismark/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        bismark_indexes_dir="resources/ref_tools/bismark/Bisulfite_Genome",
    output:
        bam="resources/ref_tools/bismark/bams/{sample}_pe_{SRA}_unsorted.bam",
        report="resources/ref_tools/bismark/bams/{sample}_{SRA}_PE_report.txt",
        nucleotide_stats="resources/ref_tools/bismark/bams/{sample}_{SRA}_pe.nucleotide_stats.txt",
    log:
        "logs/bismark/bismark_align/{sample}_{SRA}.log",
    benchmark:
        "benchmarks/Illumina_pe/bismark/bismark_align/{sample}_{SRA}.bwa.benchmark.txt"
    params:
        extra="--nucleotide_coverage",
    threads: 16
    resources:
        mem_mb=48000,
    wrapper:
        "v5.9.0/bio/bismark/bismark"


# merge bam files from different lanes
rule samtools_merge:
    input:
        get_sample_sra_bismark,
    output:
        "resources/ref_tools/bismark/bams/{sample}_pe.bam",
    log:
        "logs/bismark/samtools_merge/{sample}.log",
    params:
        extra="-n",  # optional additional parameters as string
    threads: 8
    wrapper:
        "v5.9.0/bio/samtools/merge"


rule samtools_sort:
    input:
        "resources/ref_tools/bismark/bams/{sample}_pe.bam",
    output:
        "resources/ref_tools/bismark/bams/{sample}_pe_sorted.bam",
    log:
        "logs/bismark/samtools_sort/{sample}.log",
    params:
        extra="-m 4G -n",
    threads: 8
    resources:
        mem_mb=32000,
    wrapper:
        "v5.9.0/bio/samtools/sort"


rule deduplicate_bismark:
    input:
        "resources/ref_tools/bismark/bams/{sample}_pe_sorted.bam",
    output:
        bam="resources/ref_tools/bismark/dedup/{sample}.deduplicated.bam",
        report="resources/ref_tools/bismark/dedup/{sample}.deduplication_report.txt",
    log:
        "logs/bismark/deduplicate_bismark/{sample}.log",
    params:
        extra="",  # optional params string
    benchmark:
        "benchmarks/Illumina_pe/bismark/deduplicate_bismark/{sample}.bwa.benchmark.txt"
    resources:
        mem_mb=16000,
    wrapper:
        "v5.9.0/bio/bismark/deduplicate_bismark"


# rule bismark_methylation_extractor:
#     input:
#         bam="resources/ref_tools/bismark/dedup/{sample}.deduplicated.bam",
#     output:
#         cov_zero_based="resources/ref_tools/bismark/meth/{sample}.deduplicated.bedGraph.gz.bismark.zero.cov",
#         mbias_r1="resources/ref_tools/bismark/qc/meth/{sample}.deduplicated.M-bias_R1.png",
#         # Only for PE BAMS:
#         mbias_r2="resources/ref_tools/bismark/qc/meth/{sample}.deduplicated.M-bias_R2.png",
#         mbias_report="resources/ref_tools/bismark/report/meth/{sample}.deduplicated.M-bias.txt",
#         splitting_report="resources/ref_tools/bismark/report/meth/{sample}.deduplicated_splitting_report.txt",
#         # 1-based start, 1-based end ('inclusive') methylation info: % and counts
#         methylome_CpG_cov="resources/ref_tools/bismark/meth/cov/{sample}.deduplicated.bismark.cov.gz",
#         # BedGraph with methylation percentage: 0-based start, end exclusive
#         methylome_CpG_mlevel_bedGraph="resources/ref_tools/bismark/meth/bedgraph/{sample}.deduplicated.bedGraph.gz",
#         # Primary output files: methylation status at each read cytosine position: (extremely large)
#         read_base_meth_state_cpg="resources/ref_tools/bismark/meth/CpG_context_{sample}.deduplicated.txt.gz",
#         # * You could merge CHG, CHH using: --merge_non_CpG
#         read_base_meth_state_chg="resources/ref_tools/bismark/meth/CHG_context_{sample}.deduplicated.txt.gz",
#         read_base_meth_state_chh="resources/ref_tools/bismark/meth/CHH_context_{sample}.deduplicated.txt.gz",
#         # cytosine_report="resources/ref_tools/bismark/report/meth/{sample}.deduplicated.cytosine_report.txt",
#     log:
#         "logs/bismark/bismark_methylation_extractor/{sample}.log",
#     params:
#         output_dir="resources/ref_tools/bismark/meth",  # optional output dir
#         extra="--gzip --comprehensive --bedGraph --zero_based",  # optional params string
#     benchmark:
#         "benchmarks/Illumina_pe/bismark/bismark_methylation_extractor/{sample}.bwa.benchmark.txt"
#     resources:
#         mem_mb=16000,
#     wrapper:
#         "v5.9.0/bio/bismark/bismark_methylation_extractor"


rule bismark_extract:
    input:
        bam="resources/ref_tools/bismark/dedup/{sample}.deduplicated.bam",
    output:
        cov_zero_based="resources/ref_tools/bismark/meth/{sample}.deduplicated.bedGraph.gz.bismark.zero.cov",
        # "results/single_sample/Illumina_pe/called/{sample}/result_files/CpG_context_alignment_bismark_sorted.txt",
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark/{sample}/extract_results.log",
    benchmark:
        "benchmarks/Illumina_pe/bismark/bismark_methylation_extractor/{sample}.bwa.benchmark.txt"
    resources:
        mem_mb=16000,
    shell:
        """
        mkdir -p $(dirname {output}) 2> {log}
        bismark_methylation_extractor {input} -o $(dirname {output}) --comprehensive --gzip --comprehensive --bedGraph --zero_based 2> {log}
        """


# We need this rule since the --comprehensive option in bismark_methylation_extractor
# does not create the desired bedGraph file with merged positions for forward and reverse read. We merge them manually by comparing to our candidates.
rule bismark_merge_positions:
    input:
        bedgraph="resources/ref_tools/bismark/meth/{sample}.deduplicated.bedGraph.gz.bismark.zero.cov",
        candidates=expand(
            "resources/{chrom}/candidates.bcf",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        candidates_index=expand(
            "resources/{chrom}/candidates.bcf.csi",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
    output:
        "results/single_sample/Illumina_pe/called/{sample}/result_files/bismark.bed",
    log:
        "logs/bismark/bismark_merge_positions/{sample}.log",
    conda:
        "../envs/pysam.yaml"
    benchmark:
        "benchmarks/Illumina_pe/bismark/bismark_merge_positions/{sample}.bwa.benchmark.txt"
    script:
        "../scripts/bismark_merge_positions.py"
