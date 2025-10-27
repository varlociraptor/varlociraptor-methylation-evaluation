


rule bismark_copy_chromosome:
    input:
        "resources/chromosome_{chrom}.fasta",
    output:
        "resources/ref_tools/bismark/chromosome_{chrom}.fasta",
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark/copy_chromosome_{chrom}.log",
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
        # directory("resources/ref_tools/bismark/{chrom}/Bisulfite_Genome"),
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark/prepare_genome.log",
    shell:
        """
        bismark_genome_preparation $(dirname {input}) --verbose  2> {log}
        """


rule bismark_align:
    input:
        bisulfite_folder="resources/ref_tools/bismark/Bisulfite_Genome",
        # We need this iput only for cluster execution, because else it does not find the fasta file 
        fasta=expand(
            "resources/ref_tools/bismark/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        reads1="resources/Illumina_pe/{sample, [^(?!simulated_data$)]}/{SRA}/{SRA}_1_trimmed.fastq",
        reads2="resources/Illumina_pe/{sample, [^(?!simulated_data$)]}/{SRA}/{SRA}_2_trimmed.fastq",
    output:
        "resources/ref_tools/bismark/alignment/{sample}/{SRA}/{SRA}_1_trimmed_bismark_bt2_pe.bam",
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark/{sample}/align_{SRA}.log",
    resources:
        mem_mb=32000,
    threads: 6
    shell:
        """
        mkdir -p $(dirname {output}) 2> {log}
        bismark --genome $(dirname {input.bisulfite_folder}) -1 {input.reads1} -2 {input.reads2} -o $(dirname {output}) --parallel {threads} 2> {log}
        """


# TODO: Do not accept sample=simulated_data
rule bismark_merge_bams:
    input:
        get_sample_sra_bismark,
    output:
        "resources/ref_tools/bismark/alignment/{sample}/alignment.bam",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bismark/{sample}/merge_bams.log",
    # wildcard_constraints:
    #     sample="simulated_data",
    shell:
        """
        echo {input} 2> {log}
        samtools merge -n {output} {input} 2> {log}
        """


# rule bismark_copy_chromosome_simulated:
#     input:
#         "resources/chromosome_{chrom}.fasta",
#     output:
#         "resources/ref_tools/bismark/{chrom}/chromosome_{chrom}.fasta",
#     conda:
#         "../envs/bismark.yaml"
#     log:
#         "logs/bismark/copy_chromosome_simulated_{chrom}.log",
#     shell:
#         """
#         mkdir -p $(dirname {output}) 2> {log}
#         cp {input} {output} 2> {log}
#         """


# # TODO: Gives missing output exception
# rule bismark_prepare_genome_simulated:
#     input:
#         "resources/ref_tools/bismark/{chrom}/chromosome_{chrom}.fasta",
#     output:
#         directory("resources/ref_tools/bismark/{chrom}/Bisulfite_Genome"),
#         # directory("resources/ref_tools/bismark/{chrom}/Bisulfite_Genome"),
#     conda:
#         "../envs/bismark.yaml"
#     log:
#         "logs/bismark/prepare_genome_simulated_{chrom}.log",
#     shell:
#         """
#         bismark_genome_preparation --verbose $(dirname {output}) 2> {log}
#         """


# rule bismark_align_simulated:
#     input:
#         bisulfite_folder=expand(
#             "resources/ref_tools/bismark/{chrom}/Bisulfite_Genome",
#             chrom=config["seq_platforms"].get("Illumina_pe"),
#         ),
#         chrom=expand(
#             "resources/ref_tools/bismark/{chrom}/",
#             chrom=config["seq_platforms"].get("Illumina_pe"),
#         ),
#         chromosome=expand(
#             "resources/ref_tools/bismark/{chrom}/chromosome_{chrom}.fasta",
#             chrom=config["seq_platforms"].get("Illumina_pe"),
#         ),
#         reads1=expand(
#             "resources/Illumina_pe/simulated_data/chromosome_{chrom}_f1_trimmed.fastq",
#             chrom=config["seq_platforms"].get("Illumina_pe"),
#         ),
#         reads2=expand(
#             "resources/Illumina_pe/simulated_data/chromosome_{chrom}_f2_trimmed.fastq",
#             chrom=config["seq_platforms"].get("Illumina_pe"),
#         ),
#     output:
#         "resources/ref_tools/bismark/alignment/simulated_data/chromosome_{chrom}_f1_bismark_bt2_pe.bam",
#     conda:
#         "../envs/bismark.yaml"
#     log:
#         "logs/bismark/align_simulated_{chrom}.log",
#     threads: 6
#     shell:
#         """
#         mkdir -p $(dirname {output}) 2> {log}
#         bismark --genome {input.chrom} -1 {input.reads1} -2 {input.reads2} -o $(dirname {output}) --parallel {threads} 2> {log}
#         """


# rule bismark_rename_simulated_alignment:
#     input:
#         expand(
#             "resources/ref_tools/bismark/alignment/simulated_data/chromosome_{chrom}_f1_bismark_bt2_pe.bam",
#             chrom=config["seq_platforms"].get("Illumina_pe"),
#         ),
#     output:
#         "resources/ref_tools/bismark/alignment/simulated_data/alignment.bam",
#     conda:
#         "../envs/bismark.yaml"
#     log:
#         "logs/bismark/rename_simulated_alignment.log",
#     shell:
#         """
#         mv {input} {output} 2> {log}
#         """


rule bismark_deduplicate:
    input:
        "resources/ref_tools/bismark/alignment/{sample}/alignment.bam",
    output:
        "resources/ref_tools/bismark/alignment/{sample}/alignment.deduplicated.bam",
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark/{sample}/deduplicate.log",
    shell:
        """
        deduplicate_bismark --bam {input} --output_dir $(dirname {output}) 2> {log}
        """


# It is necessary to sort by name
rule bismark_sort_bams:
    input:
        "resources/ref_tools/bismark/alignment/{sample}/alignment.deduplicated.bam",
    output:
        "resources/ref_tools/bismark/alignment/{sample}/alignment_bismark_sorted.bam",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bismark/{sample}/sort_bams.log",
    shell:
        """
        samtools merge {output} {input} 2> {log}
        samtools sort -n {input} -o {output} 2> {log}
        """


rule bismark_extract_results:
    input:
        "resources/ref_tools/bismark/alignment/{sample}/alignment_bismark_sorted.bam",
    output:
        "results/Illumina_pe/{sample}/result_files/CpG_context_alignment_bismark_sorted.txt",
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark/{sample}/extract_results.log",
    shell:
        """
        mkdir -p $(dirname {output}) 2> {log}
        bismark_methylation_extractor {input} -o $(dirname {output}) --comprehensive --merge_non_CpG 2> {log}
        """


rule bismark_to_bedGraph:
    input:
        "results/Illumina_pe/{sample}/result_files/CpG_context_alignment_bismark_sorted.txt",
    output:
        bedGraph="results/Illumina_pe/{sample}/result_files/CpG.bedGraph.gz",
        cov="results/Illumina_pe/{sample}/result_files/CpG.bismark.cov.gz",
    log:
        "logs/bismark/{sample}/to_bedGraph.log",
    wrapper:
        "v5.5.0/bio/bismark/bismark2bedGraph"


rule bismark_unzip_results:
    input:
        "results/Illumina_pe/{sample}/result_files/CpG.bismark.cov.gz",
    output:
        "results/Illumina_pe/{sample}/result_files/CpG.bismark.cov",
    conda:
        "../envs/bismark.yaml"
    log:
        "logs/bismark/{sample}/unzip_results.log",
    shell:
        """
        gunzip {input} 2> {log}
        """


rule bismark_focus_result_on_chromosome:
    input:
        "results/Illumina_pe/{sample}/result_files/CpG.bismark.cov",
    output:
        "results/single_sample/Illumina_pe/called/{sample}/result_files/bismark.bed",
    params:
        chromosome=lambda wildcards: chromosome_by_seq_platform["Illumina_pe"],
    log:
        "logs/bismark/{sample}/focus_result_on_chromosome.log",
    shell:
        """
        awk '$1 == "{params.chromosome}"' {input} > {output} 2> {log}
        """
