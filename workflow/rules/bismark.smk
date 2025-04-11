
rule copy_chromosome:
    input:
        "resources/chromosome_{chrom}.fasta",
    output:
        "resources/ref_tools/bismark/{chrom}/chromosome_{chrom}.fasta",
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        cp {input} {output}
        """


# TODO: Gives missing output exception
rule bismark_prepare:
    input:
        "resources/ref_tools/bismark/{chrom}/chromosome_{chrom}.fasta",
    output:
        directory("resources/ref_tools/bismark/{chrom}/Bisulfite_Genome"),
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        bismark_genome_preparation --verbose $(dirname {output})
        """


rule bismark_align:
    input:
        bisulfite_folder=expand(
            "resources/ref_tools/bismark/{chrom}/Bisulfite_Genome",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        chrom=expand(
            "resources/ref_tools/bismark/{chrom}/",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        chromosome=expand(
            "resources/ref_tools/bismark/{chrom}/chromosome_{chrom}.fasta",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        reads1="resources/Illumina_pe/{protocol, [^(?!simulated_data$)]}/{SRA}/{SRA}_1_trimmed.fastq",
        reads2="resources/Illumina_pe/{protocol, [^(?!simulated_data$)]}/{SRA}/{SRA}_2_trimmed.fastq",
    output:
        "resources/ref_tools/bismark/alignment/{protocol}/{SRA}/{SRA}_1_trimmed_bismark_bt2_pe.bam",
    conda:
        "../envs/bismark.yaml"
    threads: 6
    shell:
        """
        mkdir -p $(dirname {output})
        bismark --genome {input.chrom} -1 {input.reads1} -2 {input.reads2} -o $(dirname {output}) --parallel {threads}
        """


rule merge_bismark_bams:
    input:
        get_protocol_sra_bismark,
    output:
        "resources/ref_tools/bismark/alignment/{protocol}/alignment.bam",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        echo {input}
        samtools merge -n {output} {input}
        """


rule bismark_align_simulated:
    input:
        bisulfite_folder=expand(
            "resources/ref_tools/bismark/{chrom}/Bisulfite_Genome",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        chrom=expand(
            "resources/ref_tools/bismark/{chrom}/",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        chromosome=expand(
            "resources/ref_tools/bismark/{chrom}/chromosome_{chrom}.fasta",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        reads1=expand(
            "resources/Illumina_pe/simulated_data/chromosome_{chrom}_f1.fastq",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        reads2=expand(
            "resources/Illumina_pe/simulated_data/chromosome_{chrom}_f2.fastq",
            chrom=config["platforms"]["Illumina_pe"],
        ),
    output:
        "resources/ref_tools/bismark/alignment/simulated_data/chromosome_{chrom}_f1_bismark_bt2_pe.bam",
    conda:
        "../envs/bismark.yaml"
    threads: 6
    shell:
        """
        mkdir -p $(dirname {output})
        bismark --genome {input.chrom} -1 {input.reads1} -2 {input.reads2} -o $(dirname {output}) --parallel {threads}
        """


rule rename_simulated_alignment:
    input:
        expand(
            "resources/ref_tools/bismark/alignment/simulated_data/chromosome_{chrom}_f1_bismark_bt2_pe.bam",
            chrom=config["platforms"]["Illumina_pe"],
        ),
    output:
        "resources/ref_tools/bismark/alignment/simulated_data/alignment.bam",
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        mv {input} {output}
        """


rule bismark_deduplicate:
    input:
        "resources/ref_tools/bismark/alignment/{protocol}/alignment.bam",
    output:
        "resources/ref_tools/bismark/alignment/{protocol}/alignment.deduplicated.bam",
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        deduplicate_bismark --bam {input} --output_dir $(dirname {output})
        """


# It is necessary to sort by name
rule sort_bismark_bams:
    input:
        "resources/ref_tools/bismark/alignment/{protocol}/alignment.deduplicated.bam",
    output:
        "resources/ref_tools/bismark/alignment/{protocol}/alignment_bismark_sorted.bam",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        echo {input}
        samtools merge {output} {input}
        samtools sort -n {input} -o {output}
        """


rule bismark_extract_results:
    input:
        "resources/ref_tools/bismark/alignment/{protocol}/alignment_bismark_sorted.bam",
    output:
        "results/Illumina_pe/{protocol}/result_files/CpG_context_alignment_bismark_sorted.txt",
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        bismark_methylation_extractor {input} -o $(dirname {output}) --comprehensive --merge_non_CpG
        """


# rule bismark_methylation_extractor:
#     input:
#         "resources/ref_tools/bismark/alignment/{protocol}/alignment_bismark_sorted.bam",
#     output:
#         mbias_r1="qc/meth/{protocol}.M-bias_R1.png",
#         mbias_report="resources/ref_tools/bismark/meth/{protocol}.M-bias.txt",
#         splitting_report="resources/ref_tools/bismark/meth/{protocol}_splitting_report.txt",
#         # 1-based start, 1-based end ('inclusive') methylation info: % and counts
#         methylome_CpG_cov="resources/ref_tools/bismark/meth_cpg/{protocol}.bismark.cov.gz",
#         # BedGraph with methylation percentage: 0-based start, end exclusive
#         methylome_CpG_mlevel_bedGraph="resources/ref_tools/bismark/meth_cpg/{protocol}.bedGraph.gz",
#         # Primary output files: methylation status at each read cytosine position: (extremely large)
#         read_base_meth_state_cpg="resources/ref_tools/bismark/alignment/{protocol}/CpG_context.txt.gz",
#         # * You could merge CHG, CHH using: --merge_non_CpG
#         read_base_meth_state_chg="resources/ref_tools/bismark/meth/CHG_context_{protocol}.txt.gz",
#         read_base_meth_state_chh="resources/ref_tools/bismark/meth/CHH_context_{protocol}.txt.gz",
#     log:
#         "logs/meth/{protocol}.log",
#     params:
#         output_dir="meth",  # optional output dir
#         extra="--gzip --comprehensive --bedGraph",  # optional params string
#     wrapper:
#         "v5.5.0/bio/bismark/bismark_methylation_extractor"


# Example for CpG only coverage
rule bismark2bedGraph_cpg:
    input:
        "results/Illumina_pe/{protocol}/result_files/CpG_context_alignment_bismark_sorted.txt",
        # "results/Illumina_pe/{protocol}/result_files/CpG_context_test_dataset.fastq_bismark.txt",
        # "resources/ref_tools/bismark/alignment/{protocol}/CpG_context.txt.gz",
    output:
        bedGraph="results/Illumina_pe/{protocol}/result_files/CpG.bedGraph.gz",
        cov="results/Illumina_pe/{protocol}/result_files/CpG.bismark.cov.gz",
    log:
        "logs/meth_cpg/{protocol}_CpG.log",
    wrapper:
        "v5.5.0/bio/bismark/bismark2bedGraph"


rule unzip_bismark_results:
    input:
        "results/Illumina_pe/{protocol}/result_files/CpG.bismark.cov.gz",
    output:
        # "results/Illumina_pe/{protocol}/result_files/alignment_bismark_sorted.bedGraph",
        "results/Illumina_pe/{protocol}/result_files/CpG.bismark.cov",
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        gunzip {input}
        """


rule bismark_focus_on_chromosome:
    input:
        "results/Illumina_pe/{protocol}/result_files/CpG.bismark.cov",
    output:
        "results/Illumina_pe/{protocol}/result_files/bismark.bed",
    params:
        chromosome=lambda wildcards: chromosome_by_platform["Illumina_pe"],
    shell:
        """
        awk '$1 == "{params.chromosome}"' {input} > {output}
        """
