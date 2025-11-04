# https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/


rule bismark_copy_genome:
    input:
        "resources/genome.fasta",
    output:
        "resources/ref_tools/bismark/genome.fasta",
    log:
        "logs/bismark/copy_chromosome.log",
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
        "logs/bismark/copy_chromosome_{chrom}.log",
    shell:
        """
        mkdir -p $(dirname {output}) 2> {log}
        cp {input} {output} 2> {log}
        """


# # Prepare genome and indices
# rule bismark_genome_preparation_fa:
#     input:
#         expand(
#             "resources/ref_tools/bismark/genome.fasta",
#             chrom=config["seq_platforms"].get("Illumina_pe"),
#         ),
#     #         "resources/genome.fasta",
#     output:
#         directory("resources/ref_tools/bismark/Bisulfite_Genome"),
#     log:
#         "logs/Bisulfite_Genome.log",
#     params:
#         extra="",  # optional params string
#     wrapper:
#         "v5.9.0/bio/bismark/bismark_genome_preparation"


# # TODO: Gives missing output exception
# rule bismark_prepare_genome:
#     input:
#         "resources/ref_tools/bismark/genome.fasta",
#     # "resources/genome.fasta",
#     output:
#         directory("resources/ref_tools/bismark/Bisulfite_Genome"),
#         # directory("resources/ref_tools/bismark/{chrom}/Bisulfite_Genome"),
#     conda:
#         "../envs/bismark.yaml"
#     log:
#         "logs/bismark/prepare_genome.log",
#     shell:
#         """
#         bismark_genome_preparation $(dirname {input}) --verbose  2> {log}
#         """


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
        "logs/bams/{sample}_{SRA}.log",
    benchmark:
        "benchmarks/Illumina_pe/bismark/bismark_align/{sample}_{SRA}.bwa.benchmark.txt"
    params:
        extra="--nucleotide_coverage",
    threads: 16
    resources:
        mem_mb=48000  
    wrapper:
        "v5.9.0/bio/bismark/bismark"


# merge bam files from different lanes
rule samtools_merge:
    input:
        get_sample_sra_bismark,
    output:
        "resources/ref_tools/bismark/bams/{sample}_pe.bam",
    log:
        "logs/bams/{sample}.log",
    params:
        extra="-n",  # optional additional parameters as string
    threads: 8
    wrapper:
        "v5.9.0/bio/samtools/merge"


rule samtools_sort:
    input:
        "resources/ref_tools/bismark/bams/{sample}_pe.bam",

        # "resources/ref_tools/bismark/dedup/{sample}.deduplicated.bam",
    output:
        "resources/ref_tools/bismark/bams/{sample}_pe_sorted.bam",
        # "resources/ref_tools/bismark/bams/{sample}_pe_name_sorted.bam",
    log:
        "logs/bams/{sample}_pe_name_sorted.log",
    params:
        extra="-m 4G -n",
    threads: 8
    resources:
        mem_mb=32000   # ✅ 32GB
    wrapper:
        "v5.9.0/bio/samtools/sort"


rule deduplicate_bismark:
    input:
        "resources/ref_tools/bismark/bams/{sample}_pe_sorted.bam",
    output:
        bam="resources/ref_tools/bismark/dedup/{sample}.deduplicated.bam",
        report="resources/ref_tools/bismark/dedup/{sample}.deduplication_report.txt",
    log:
        "logs/dedup/{sample}.deduplicated.log",
    params:
        extra="",  # optional params string
    benchmark:
        "benchmarks/Illumina_pe/bismark/deduplicate_bismark/{sample}.bwa.benchmark.txt"
    resources:
        mem_mb=16000  
    wrapper:
        "v5.9.0/bio/bismark/deduplicate_bismark"





rule bismark_methylation_extractor:
    input:
        bam="resources/ref_tools/bismark/dedup/{sample}.deduplicated.bam",

        # "resources/ref_tools/bismark/bams/{sample}_pe_name_sorted.bam",
    output:
        cov_zero_based="resources/ref_tools/bismark/meth/{sample}.deduplicated.bedGraph.gz.bismark.zero.cov",
        mbias_r1="resources/ref_tools/bismark/qc/meth/{sample}.deduplicated.M-bias_R1.png",
        # Only for PE BAMS:
        mbias_r2="resources/ref_tools/bismark/qc/meth/{sample}.deduplicated.M-bias_R2.png",
        mbias_report="resources/ref_tools/bismark/report/meth/{sample}.deduplicated.M-bias.txt",
        splitting_report="resources/ref_tools/bismark/report/meth/{sample}.deduplicated_splitting_report.txt",
        # 1-based start, 1-based end ('inclusive') methylation info: % and counts
        methylome_CpG_cov="resources/ref_tools/bismark/meth/cov/{sample}.deduplicated.bismark.cov.gz",
        # BedGraph with methylation percentage: 0-based start, end exclusive
        methylome_CpG_mlevel_bedGraph="resources/ref_tools/bismark/meth/bedgraph/{sample}.deduplicated.bedGraph.gz",
        # Primary output files: methylation status at each read cytosine position: (extremely large)
        read_base_meth_state_cpg="resources/ref_tools/bismark/meth/CpG_context_{sample}.deduplicated.txt.gz",
        # * You could merge CHG, CHH using: --merge_non_CpG
        read_base_meth_state_chg="resources/ref_tools/bismark/meth/CHG_context_{sample}.deduplicated.txt.gz",
        read_base_meth_state_chh="resources/ref_tools/bismark/meth/CHH_context_{sample}.deduplicated.txt.gz",
        # cytosine_report="resources/ref_tools/bismark/report/meth/{sample}.deduplicated.cytosine_report.txt",
    log:
        "logs/meth/{sample}.log",
    params:
        output_dir="resources/ref_tools/bismark/meth",  # optional output dir
        extra="--gzip --comprehensive --bedGraph --zero_based",  # optional params string
    benchmark:
        "benchmarks/Illumina_pe/bismark/bismark_methylation_extractor/{sample}.bwa.benchmark.txt"
    resources:
        mem_mb=16000
    wrapper:
        "v5.9.0/bio/bismark/bismark_methylation_extractor"


# # Example for CpG only coverage
# rule bismark2bedGraph_cpg:
#     input:
#         "resources/ref_tools/bismark/meth/CpG_context_{sample}.txt.gz",
#     output:
#         bedGraph="meth_cpg/{sample}_CpG.bedGraph.gz",
#         cov="meth_cpg/{sample}_CpG.bismark.cov.gz",
#     log:
#         "logs/meth_cpg/{sample}_CpG.log",
#     wrapper:
#         "v7.3.0/bio/bismark/bismark2bedGraph"


# We need this rule since the --comprehensive option in bismark_methylation_extractor
# does not create the desired bedGraph file with merged positions for forward and reverse read. We merge them manually by comparing to our candidates.
rule bismark_merge_positions:
    input:
        bedgraph="resources/ref_tools/bismark/meth/{sample}.deduplicated.bedGraph.gz.bismark.zero.cov",
        candidates=expand(
            "resources/{chrom}/candidates.vcf",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
    output:
        "results/single_sample/Illumina_pe/called/{sample}/result_files/bismark.bed",
    log:
        "logs/bismark/{sample}/bismark_merge_positions/bismark.log",
    benchmark:
        "benchmarks/Illumina_pe/bismark/bismark_merge_positions/{sample}.bwa.benchmark.txt"
    script:
        "../scripts/bismark_merge_positions.py"


# rule bismark_unzip:
# input:
#     "resources/ref_tools/bismark/meth/bedgraph/{sample}.bedGraph.gz",
# output:
#     "results/single_sample/Illumina_pe/called/{sample}/result_files/bismark.bed",
# log:
#     "logs/bismark/{sample}/unzip_resources/bismark.log",
# shell:
#     """
#     mkdir -p $(dirname {output})
#     gunzip -c {input} > {output} 2> {log}
#     """
# rule bismark_focus_result_on_chromosome:
#     input:
#         "resources/ref_tools/bismark/single_sample/Illumina_pe/called/{sample}/result_files/CpG.bismark.cov",
#     output:
#         "resources/ref_tools/bismark/single_sample/Illumina_pe/called/{sample}/result_files/bismark.bed",
#     params:
#         chromosome=lambda wildcards: chromosome_by_seq_platform["Illumina_pe"],
#     log:
#         "logs/bismark/{sample}/focus_result_on_chromosome.log",
#     shell:
#         """
#         awk '$1 == "{params.chromosome}"' {input} > {output} 2> {log}
#         """
# rule merge_cpg_symmetrically:
#     input:
#         cpg="resources/ref_tools/bismark/meth/CpG_context_{sample}.txt.gz",
#     output:
#         "resources/ref_tools/bismark/meth/CpG_context_merged/{sample}_CpG_merged.txt.gz",
#     log:
#         "logs/meth/{sample}_CpG_merged.log",
#     conda:
#         "../envs/bismark.yaml"
#     shell:
#         """
#         bismark2cytosine \
#             --merge_CpG \
#             --gzip \
#             --output {output} \
#             {input.cpg} \
#             > {log} 2>&1
#         """
