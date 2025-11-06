ref_gene = config.get("sample", {})
chromosomes = set(chromosome for chromosome in config["seq_platforms"].values())


rule download_genome:
    output:
        "resources/genome.fasta",
    params:
        species=ref_gene.get("species"),
        datatype=ref_gene.get("datatype"),
        build=ref_gene.get("build"),
        release=ref_gene.get("release"),
    log:
        "logs/download_data/download_genome/download.log",
    cache: "omit-software"
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule genome_index:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/download_data/genome_index/index.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input} 2> {log}"


rule focus_genome_on_chromosome:
    input:
        "resources/genome.fasta",
    output:
        "resources/chromosome_{chromosome}.fasta",
    log:
        "logs/download_data/focus_genome_on_chromosome/{chromosome}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    shell:
        "samtools faidx {input} {wildcards.chromosome} > {output} 2> {log}"


rule chromosome_index:
    input:
        "resources/chromosome_{chromosome}.fasta",
    output:
        "resources/chromosome_{chromosome}.fasta.fai",
    log:
        "logs/download_data/chromosome_index/{chromosome}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input} 2> {log}"


rule rename_chromosome_in_fasta:
    input:
        expand(
            "resources/chromosome_{chromosome}.fasta",
            chromosome=[chr for chr in chromosomes],
        ),
    output:
        "resources/chr_chromosome_{chromosome}.fasta",
    log:
        "logs/download_data/rename_chromosome_in_fasta/{chromosome}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_chrom_in_fasta.py"


rule get_fastq_pe:
    output:
        "resources/Illumina_pe/{sample,(?!simulated_data).*}/{SRA}/{SRA}_1.fastq",
        "resources/Illumina_pe/{sample,(?!simulated_data).*}/{SRA}/{SRA}_2.fastq",
    log:
        "logs/download_data/get_fastq_pe/{sample}_{SRA}_{SRA}.log",
    params:
        extra="--skip-technical",
    threads: 6
    # wildcard_constraints:
    #     sample="",
    wrapper:
        "v7.1.0/bio/sra-tools/fasterq-dump"


rule get_fastq_se:
    output:
        "resources/Illumina_se/{sample}/{SRA}/{SRA}.fastq",
    log:
        "logs/download_data/get_fastq_se/{sample}_{SRA}_{SRA}.log",
    params:
        extra="--skip-technical",
    threads: 6
    wrapper:
        "v7.1.0/bio/sra-tools/fasterq-dump"


rule trim_fastq_pe:
    input:
        first="resources/Illumina_pe/{sample,(?!simulated_data).*}/{SRA}/{SRA}_1.fastq",
        second="resources/Illumina_pe/{sample,(?!simulated_data).*}/{SRA}/{SRA}_2.fastq",
    output:
        first="resources/Illumina_pe/{sample,(?!simulated_data).*}/{SRA}/{SRA}_1_trimmed.fastq",
        second="resources/Illumina_pe/{sample,(?!simulated_data).*}/{SRA}/{SRA}_2_trimmed.fastq",
    log:
        "logs/download_data/trim_fastq_pe/{sample}_{SRA}_{SRA}.log",
    conda:
        "../envs/fastp.yaml"
    wildcard_constraints:
        sample="^(?!simulated_data$).*",
    shell:
        "fastp --in1 {input.first} --in2 {input.second} --out1 {output.first} --out2 {output.second} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis 2> {log}"


rule trim_fastq_se:
    input:
        first="resources/Illumina_se/{sample}/{SRA}/{SRA}.fastq",
    output:
        first="resources/Illumina_se/{sample}/{SRA}/{SRA}_trimmed.fastq",
    log:
        "logs/download_data/trim_fastq_se/{sample}_{SRA}_{SRA}.log",
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp --in1 {input.first} --out1 {output.first} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis 2> {log}"


rule get_pacbio_data:
    output:
        alignment="resources/PacBio/{sample}/{SRA}/alignment.bam",
    params:
        url=lambda wildcards: config.get(str(wildcards.SRA)),
        chromosome=f"chr{config['seq_platforms'].get('PacBio')}",
    log:
        "logs/download_data/get_pacbio_data/{sample}_{SRA}.log",
    resources:
        mem_mb=4096,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -b {params.url} {params.chromosome} > {output.alignment} 2> {log}"


# TODO: Does not work for replicate2. You have to download this manually with wget right now
rule get_nanopore_data:
    output:
        alignment="resources/Nanopore/{sample}/{SRA}/alignment.bam",
    params:
        url=lambda wildcards: config.get(str(wildcards.SRA)),
        chromosome=f"chr{config['seq_platforms'].get('Nanopore')}",
    log:
        "logs/download_data/get_nanopore_data/{sample}_{SRA}.log",
    resources:
        mem_mb=4096,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -b {params.url} {params.chromosome} > {output.alignment} 2> {log}"
