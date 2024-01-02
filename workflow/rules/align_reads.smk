rule index_genome:
    input:
        "resources/genome.fasta"
    output:
        "resources/genome.fasta.bwameth.c2t"
    conda:
        "../envs/bwa-meth.yaml"
    shell:
        "bwameth.py index-mem2 {input}"

rule align_reads:
    input:
        fasta_index="resources/genome.fasta.bwameth.c2t",
        fasta="resources/genome.fasta",
        reads1="resources/Illumina/{protocol}/{SRA}/{SRA}_1_trimmed.fastq",
        reads2="resources/Illumina/{protocol}/{SRA}/{SRA}_2_trimmed.fastq",
    output:
        "resources/Illumina/{protocol}/{SRA}/alignment.bam"
    conda:
        "../envs/bwa-meth.yaml"
    log:
        # "logs/align_reads_{scatteritem}{SRA}.log",
        "logs/align_reads{SRA}Illumina/{protocol}.log",
    threads: 30
    resources:
        mem_mb=512,
    shell:
        """
        bwameth.py --threads {threads} --reference {input.fasta} {input.reads1} {input.reads2}  | samtools view -S -b - > {output}
        """

rule sort_aligned_reads:
    input:
        "resources/{platform}/{protocol}/{SRA}/alignment.bam"
    output:
        "resources/{platform}/{protocol}/{SRA}/alignment_sorted.bam"
    log:
        "logs/sort_aligned_reads{SRA}{platform}/{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        samtools sort -@ {threads}  {input} -o {output}    
        """

rule aligned_reads_index:
    input:
        "resources/{platform}/{protocol}/{SRA}/alignment_sorted.bam"
    output:
        "resources/{platform}/{protocol}/{SRA}/alignment_sorted.bam.bai",
    log:
        "logs/aligned_reads_to_bam{platform}/{protocol}{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        samtools index -@ {threads} {params.pipeline_path}{input}
      """


# params.chromosome muss mak 21 und mal chr21 sein (Illumina 21, PacBio chr21)
rule focus_aligned_reads_chrom:
    input:
        bam="resources/{platform}/{protocol}/{SRA}/alignment_sorted.bam",
        index="resources/{platform}/{protocol}/{SRA}/alignment_sorted.bam.bai",
    output:
        bam="resources/{platform}/{protocol}/{SRA}/alignment_focused.bam"
    log:
        "logs/chromosome_index{SRA}{platform}/{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        chromosome=lambda wildcards: "chr" + chromosome_conf["chromosome"] if wildcards.platform == "PacBio" or wildcards.platform == "Nanopore" else chromosome_conf["chromosome"]
    threads: 10
    shell:
        """ 
        samtools view -b -o {output.bam} {input} {params.chromosome}
        """

rule filter_mapping_quality:
    input:
        "resources/{platform}/{protocol}/{SRA}/alignment_focused.bam"
    output:
        "resources/{platform}/{protocol}/{SRA}/alignment_focused_filtered.bam"
    log:
        "logs/filter_mappingq{platform}/{protocol}{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        min_quality=10
    threads: 10
    shell:
        """
        samtools view -q {params.min_quality} -b -o {output} {input}
        """



rule markduplicates_bam:
    input:
        bams="resources/{platform}/{protocol}/{SRA}/alignment_focused_filtered.bam"
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="resources/{platform}/{protocol}/{SRA}/alignment_focused_dedup.bam",
        metrics="resources/{platform}/{protocol}/{SRA}/alignment_focused_dedup.metrics.txt"
    log:
        "logs/dedup_bam/{SRA}{platform}/{protocol}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v2.6.0/bio/picard/markduplicates"


# rule aligned_reads_dedup_index:
#     input:
#         bam="resources/{platform}/{protocol}/{SRA}/alignment_focused_dedup.bam",
#     output:
#         bam="resources/{platform}/{protocol}/{SRA}/alignment_focused_dedup.bam.bai",
#     log:
#         "logs/aligned_reads_to_bam{platform}/{protocol}{SRA}.log",
#     conda:
#         "../envs/samtools.yaml"
#     params:
#         pipeline_path=config["pipeline_path"],
#     threads: 10
#     shell:
#         """
#         samtools index -@ {threads} {params.pipeline_path}{input}
#         """



def get_protocol_sra(wildcards):
    
    platform = wildcards.platform
    protocol = wildcards.protocol
    accession_numbers = config["data"][platform][protocol]
    return ["resources/" + platform + "/" + protocol + "/" + SRA + "/alignment_focused_dedup.bam" for SRA in accession_numbers]



rule merge_bams:
    input:
        get_protocol_sra
    output:
        "resources/{platform, [^/]+}/{protocol, [^/]+}/alignment_focused_dedup.bam"
    log:
        "logs/aligned_reads_to_bam{platform}/{protocol}.log",
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

rule downsample_bams:
    input:
        "resources/{platform}/{protocol}/alignment_focused_dedup.bam"
    output:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup.bam"
    log:
        "logs/downsample_bams{platform}/{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        samtools view -s 0.99 -b -o {output} {input}
        """

rule bam_to_sam:
    input:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup.bam"
    output:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup.sam"
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads:
        10
    shell:
        """
        samtools view -@ {threads} -h  -o {params.pipeline_path}/{output}  {params.pipeline_path}{input}   

        """


rule rename_chromosomes_in_sam:
    input:
        bam="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup.sam"
    output:
        renamed="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.sam"
    script:
        "../scripts/rename_chromosomes.py"

# Funktioniert nicht mit snakemake aber manuell...
rule sam_to_bam:
    input:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.sam"
    output:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam"
    conda:
        "../envs/samtools.yaml"
    threads:
        10
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """
        samtools view -bS {params.pipeline_path}/{input} > {params.pipeline_path}{output}   
        """


rule aligned_downsampled_reads_dedup_index:
    input:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam"
    output:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai"
    log:
        "logs/aligned_reads_to_bam{platform}/{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        samtools index -@ {threads} {params.pipeline_path}{input}
        """