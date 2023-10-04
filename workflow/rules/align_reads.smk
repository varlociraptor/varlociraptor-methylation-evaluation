accession = config["accession"]
read_type = config["read_type"]


aligned_bam_reads = (
    "resources/{SRA}/aligned-reads-pacbio.bam"
    if read_type == "PacBio"
    else "resources/{SRA}/aligned-reads-illumina.bam"
)

aligned_bam_reads_index = (
    "resources/{SRA}/aligned-reads-pacbio-sorted.bam.bai"
    if read_type == "PacBio"
    else "resources/{SRA}/aligned-reads-illumina-sorted.bam.bai"
)

aligned_sam_reads = (
    "resources/{SRA}/aligned-reads-pacbio.sam"
    if read_type == "PacBio"
    else "resources/{SRA}/aligned-reads-illumina.sam"
)

aligned_bam_reads_sorted = (
    "resources/{SRA}/aligned-reads-pacbio-sorted.bam"
    if read_type == "PacBio"
    else "resources/{SRA}/aligned-reads-illumina-sorted.bam"
)

aligned_sam_reads_sorted = (
    "resources/{SRA}/aligned-reads-pacbio-sorted.sam"
    if read_type == "PacBio"
    else "resources/{SRA}/aligned-reads-illumina-sorted.sam"
)



rule align_reads:
    input:
        fasta="resources/genome.fasta",
        reads1="resources/{SRA}/{SRA}_1.fastq",
        reads2="resources/{SRA}/{SRA}_2.fastq",
    output:
        aligned_bam_reads
    conda:
        "../envs/bwa-meth.yaml"
    log:
        # "logs/align_reads_{scatteritem}{SRA}.log",
        "logs/align_reads{SRA}.log",
    threads: 30
    shell:
        """
        bwameth.py index-mem2 {input.fasta}
        bwameth.py --threads {threads} --reference {input.fasta} {input.reads1} {input.reads2}  | samtools view -S -b - > {output}
        """


rule sort_aligned_reads:
    input:
        aligned_bam_reads,
    output:
        temp(aligned_bam_reads_sorted),
    log:
        "logs/sort_aligned_reads{SRA}.log",
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
        aligned_bam_reads_sorted,
    output:
        temp(aligned_bam_reads_index),
    log:
        "logs/aligned_reads_to_bam{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        samtools index -@ {threads} {params.pipeline_path}/{input}
        """

rule filter_aligned_reads:
    input:
        bam=aligned_bam_reads_sorted,
        index=aligned_bam_reads_index,
    output:
        "resources/{SRA}/alignment_focused.bam"
    log:
        "logs/chromosome_index{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        chromosome=chromosome_conf["chromosome"],
    threads: 10
    shell:
        """ 
        samtools view -b -o {output} {input.bam} {params.chromosome}
        """



rule aligned_reads_filtered_index:
    input:
        "resources/{SRA}/alignment_focused.bam"
    output:
        temp("resources/{SRA}/alignment_focused.bam.bai")
    log:
        "logs/aligned_reads_to_bam{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        samtools index -@ {threads} {params.pipeline_path}/{input}
        """


