rule candidates_to_vcf:
    input:
        "resources/candidates.bcf",
    output:
        "resources/candidates.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/convert_to_vcf.log",
    threads: 10
    shell:
        """
        bcftools view --threads {threads} {input} > {output}
        """

rule candidate_splits_to_vcf:
    input:
        "resources/candidates_{scatteritem}.bcf"
    output:
        temp("resources/candidates_{scatteritem}.vcf"),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/convert_splits_to_vcf_{scatteritem}.log",
    shell:
        """
        bcftools view {input} > {output}
        """

rule aligned_reads_sorted_sam:
    input:
        "resources/{SRA}/aligned-reads-illumina-sorted.bam",
    output:
        "resources/{SRA}/aligned-reads-illumina-sorted.sam",
    log:
        "logs/aligned_reads_sorted_sam{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """ 
        samtools view -@ {threads} -h  -o {params.pipeline_path}/{output}  {params.pipeline_path}{input}   
        """

rule aligned_reads_filtered_sam:
    input:
        "resources/{SRA}/alignment_focused.bam",
    output:
        "resources/{SRA}/alignment_focused.sam",
    log:
        "logs/aligned_reads_filtered_sam{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """ 
        samtools view -@ {threads} -h  -o {params.pipeline_path}/{output}  {params.pipeline_path}{input}   
        """


rule observations_to_vcf:
    input:
        "results/{SRA}/normal_{scatteritem}.bcf",
        # "results/{SRA}/normal.bcf",
    output:
        temp("results/{SRA}/normal_{scatteritem}.vcf"),
        # "results/{SRA}/normal.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/convert_to_vcf_{scatteritem}{SRA}.log",
    threads: 10
    shell:
        """
        bcftools view --threads {threads} {input} > {output}
        """

rule aligned_reads_complete_sam:
    input:
        "resources/{protocol}/alignment_focused_downsampled_dedup.bam",
    output:
        "resources/{protocol}/alignment_focused_downsampled_dedup.sam",
    log:
        "logs/aligned_reads_sorted_sam{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """ 
        samtools view -@ {threads} -h  -o {params.pipeline_path}/{output}  {params.pipeline_path}{input}   
        """