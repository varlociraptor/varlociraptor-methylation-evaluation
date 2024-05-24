# rule bismark:
#     input:
#         genome="resources/genome.fasta",
#         reads1="resources/Illumina_pe/{protocol}/{SRA}/{SRA}_1_trimmed.fastq",
#         reads2="resources/Illumina_pe/{protocol}/{SRA}/{SRA}_2_trimmed.fastq",
#     output:
#         "results/ref_tools/bismark/{protocol}/{SRA}_trimmed_bismark_bt2_pe.deduplicated.bedGraph",
#         "resources/ref_tools/bismark/results/{protocol}/{SRA}_trimmed_bismark_bt2_pe.deduplicated.bedGraph",
# conda:
#     "../envs/bismark.yaml"
#     params:
#         bismark_dir="~/Documents/Promotion/varlociraptor-methylation-evaluation/resources/ref_tools/bismark",
#         # bismark_dir=config["pipeline_path"] + "resources/ref_tools/bismark",
#     threads: 6
#     shell:
#         """
#         mkdir -p {params.bismark_dir}
#         cp {input.genome} {params.bismark_dir}
#         bismark_genome_preparation --verbose {params.bismark_dir}
#         mkdir -p {params.bismark_dir}/alignment/{wildcards.protocol}
#         bismark --genome {params.bismark_dir} -1 {input.reads1} -2 {input.reads2} -o {params.bismark_dir}/alignment/{wildcards.protocol} --parallel {threads}
#         deduplicate_bismark --bam {params.bismark_dir}/alignment/{wildcards.protocol}/*.bam --output_dir {params.bismark_dir}/alignment/{wildcards.protocol}
#         mkdir -p {params.bismark_dir}/results/{wildcards.protocol}
#         bismark_methylation_extractor --bedGraph {params.bismark_dir}/alignment/{wildcards.protocol}/*.deduplicated.bam -o {params.bismark_dir}/results/{wildcards.protocol}
#         gunzip {params.bismark_dir}/results/{wildcards.protocol}/*.bedGraph.gz
#         """


rule bismark_prepare:
    input:
        "resources/genome.fasta",
    output:
        directory("resources/ref_tools/bismark/Bisulfite_Genome"),
    params:
        bismark_dir=config["pipeline_path"] + "resources/ref_tools/bismark",
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        mkdir -p {params.bismark_dir}
        cp {input} {params.bismark_dir}
        bismark_genome_preparation --verbose {params.bismark_dir}
        """


rule bismark_align:
    input:
        "resources/ref_tools/bismark/Bisulfite_Genome",
        reads1="resources/Illumina_pe/{protocol}/{SRA}/{SRA}_1_trimmed.fastq",
        reads2="resources/Illumina_pe/{protocol}/{SRA}/{SRA}_2_trimmed.fastq",
    output:
        "resources/ref_tools/bismark/alignment/{protocol}/{SRA}/{SRA}_1_trimmed_bismark_bt2_pe.bam",
    params:
        bismark_dir=config["pipeline_path"] + "resources/ref_tools/bismark",
    conda:
        "../envs/bismark.yaml"
    threads: 6
    shell:
        """
        mkdir -p $(dirname {output})
        bismark --genome {params.bismark_dir} -1 {input.reads1} -2 {input.reads2} -o $(dirname {output}) --parallel {threads}
        """


rule bismark_deduplicate:
    input:
        "resources/ref_tools/bismark/alignment/{protocol, [^/]+}/{SRA}/{SRA}_1_trimmed_bismark_bt2_pe.bam",
    output:
        "resources/ref_tools/bismark/alignment/{protocol, [^/]+}/{SRA}/{SRA}_1_trimmed_bismark_bt2_pe.deduplicated.bam",
    params:
        bismark_dir=config["pipeline_path"] + "resources/ref_tools/bismark",
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        deduplicate_bismark --bam {input} --output_dir $(dirname {output})
        """


rule merge_bismark_bams:
    input:
        get_protocol_sra_bismark,
    output:
        "resources/ref_tools/bismark/alignment/{protocol}/alignment_bismark.bam",
    log:
        "logs/merge_bams_{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        echo {input}
        samtools merge {output} {input}
        """


# It is necessary to sort by name
rule sort_bismark_bams:
    input:
        "resources/ref_tools/bismark/alignment/{protocol}/alignment_bismark.bam",
    output:
        "resources/ref_tools/bismark/alignment/{protocol}/alignment_bismark_sorted.bam",
    log:
        "logs/merge_bams_{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
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
        "results/ref_tools/bismark/{protocol}/alignment_bismark_sorted.bedGraph.gz",
    params:
        bismark_dir=config["pipeline_path"] + "resources/ref_tools/bismark",
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        bismark_methylation_extractor --bedGraph {input} -o $(dirname {output}) --split_by_chromosome
        """


rule unzip_bismark_results:
    input:
        "results/ref_tools/bismark/{protocol}/alignment_bismark_sorted.bedGraph.gz",
    output:
        "results/ref_tools/bismark/{protocol}/alignment_bismark_sorted.bedGraph",
    params:
        bismark_dir=config["pipeline_path"] + "resources/ref_tools/bismark",
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        gunzip {input}
        """


rule bismark_focus_on_chromosome:
    input:
        "results/ref_tools/bismark/{protocol}/alignment_bismark_sorted.bedGraph",
    output:
        "results/ref_tools/bismark/{protocol}/bismark.bed",
    params:
        chromosome=lambda wildcards: chromosome_by_platform["Illumina_pe"],
    shell:
        """
        awk '$1 == "{params.chromosome}"' {input} > {output}
        """
