scattergather:
    split_candidates=20,


rule candidates_to_vcf:
    input:
        "resources/{chro}/candidates.bcf",
    output:
        "resources/{chro}/candidates.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/{chro}/candidates_to_vcf.log",
    threads: 10
    shell:
        """
        bcftools view --threads {threads} {input} > {output}
        """


rule candidate_splits_to_vcf:
    input:
        "resources/{chro}/candidates_{scatteritem}.bcf",
    output:
        temp("resources/{chro}/candidates_{scatteritem}.vcf"),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/{chro}/candidate_splits_to_vcf_{scatteritem}.log",
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
        "logs/aligned_reads_sorted_sam_{SRA}.log",
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
        "logs/aligned_reads_filtered_sam_{SRA}.log",
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
    output:
        temp("results/{SRA}/normal_{scatteritem}.vcf"),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/observations_to_vcf_{SRA}_{scatteritem}.log",
    threads: 10
    shell:
        """
        bcftools view --threads {threads} {input} > {output}
        """


rule bcf_indices:
    input:
        "results/{platform}/{protocol}/normal_{scatteritem}.bcf",
    output:
        "results/{platform}/{protocol}/normal_{scatteritem}.bcf.csi",
    log:
        "logs/split_candidates{platform}_{protocol}_{scatteritem}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """
        bcftools index  {params.pipeline_path}{input}
        """


rule single_bcf:
    input:
        bcf=gather.split_candidates(
            "results/{{platform}}/{{protocol}}/normal_{scatteritem}.bcf"
        ),
        indices=gather.split_candidates(
            "results/{{platform}}/{{protocol}}/normal_{scatteritem}.bcf.csi"
        ),
    output:
        "results/{platform}/{protocol}/normal.bcf",
    log:
        "logs/single_bcf_{platform}_{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """
        bcftools merge -o {params.pipeline_path}{output} {params.pipeline_path}{input.bcf}
        """


rule normal_to_vcf:
    input:
        "results/{platform}/{protocol}/normal.bcf",
    output:
        "results/{platform}/{protocol}/normal.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/normal_to_vcf_{platform}_{protocol}.log",
    threads: 10
    shell:
        """
        bcftools view --threads {threads} {input} > {output}
        """
