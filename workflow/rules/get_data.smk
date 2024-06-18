ref_gene = config["sample"]
chromosomes = set(chromosome for chromosome in config["platforms"].values())


rule get_genome:
    output:
        "resources/genome.fasta",
    params:
        species=ref_gene["species"],
        datatype=ref_gene["datatype"],
        build=ref_gene["build"],
        release=ref_gene["release"],
    log:
        "logs/get_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule get_chromosome_from_genome:
    input:
        "resources/genome.fasta",
    output:
        "resources/chromosome_{chromosome}.fasta",
    log:
        "logs/get_chromosome_from_genome_{chromosome}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    shell:
        """ 
        samtools faidx {input} {wildcards.chromosome} > {output}
        """


rule chromosome_index:
    input:
        "resources/chromosome_{chromosome}.fasta",
    output:
        "resources/chromosome_{chromosome}.fasta.fai",
    log:
        "logs/chromosome_index_{chromosome}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        samtools faidx {params.pipeline_path}{input}
        """


rule add_chr_to_fasta:
    input:
        expand(
            "resources/chromosome_{chromosome}.fasta",
            chromosome=[chr for chr in chromosomes],
        ),
    output:
        "resources/chr_chromosome_{chromosome}.fasta",
    log:
        "logs/add_chr_to_fasta_{chromosome}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_chromosome.py"


rule genome_index:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome_index.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        samtools faidx {params.pipeline_path}{input}
        """


rule get_fastq_pe:
    output:
        "resources/Illumina_pe/{protocol}/{SRA}/{accession}_1.fastq",
        "resources/Illumina_pe/{protocol}/{SRA}/{accession}_2.fastq",
    log:
        "logs/get_fastq_pe_{protocol}_{SRA}_{accession}.log",
    params:
        extra="--skip-technical",
    threads: 6
    conda:
        "../envs/fastq-wrapper.yaml"
    wrapper:
        "v3.0.2/bio/sra-tools/fasterq-dump"


rule get_fastq_se:
    output:
        "resources/Illumina_se/{protocol}/{SRA}/{accession}.fastq",
    log:
        "logs/get_fastq_se_{protocol}_{SRA}_{accession}.log",
    params:
        extra="--skip-technical",
    threads: 6
    wrapper:
        "v3.3.3/bio/sra-tools/fasterq-dump"


rule trim_fastq_pe:
    input:
        first="resources/Illumina_pe/{protocol}/{SRA}/{accession}_1.fastq",
        second="resources/Illumina_pe/{protocol}/{SRA}/{accession}_2.fastq",
    output:
        first="resources/Illumina_pe/{protocol}/{SRA}/{accession}_1_trimmed.fastq",
        second="resources/Illumina_pe/{protocol}/{SRA}/{accession}_2_trimmed.fastq",
    log:
        "logs/trim_fastq_pe_{protocol}_{SRA}_{accession}.log",
    conda:
        "../envs/fastp.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        fastp --in1 {input.first} --in2 {input.second} --out1 {output.first} --out2 {output.second} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis
        """


rule trim_fastq_se:
    input:
        first="resources/Illumina_se/{protocol}/{SRA}/{accession}.fastq",
    output:
        first="resources/Illumina_se/{protocol}/{SRA}/{accession}_trimmed.fastq",
    log:
        "logs/trim_fastq_se_{protocol}_{SRA}_{accession}.log",
    conda:
        "../envs/fastp.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        fastp --in1 {input.first} --out1 {output.first} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis
        """


rule get_pacbio_data:
    output:
        alignment="resources/PacBio/{protocol}/{SRA}/alignment.bam",
    params:
        url=lambda wildcards: config[str(wildcards.SRA)],
    log:
        "logs/get_pacbio_data_{protocol}_{SRA}.log",
    resources:
        mem_mb=4096,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_pacbio_data.py"


rule get_nanopore_header:
    output:
        header="resources/Nanopore/{protocol}/{SRA}/header.sam",
    params:
        pipeline_path=config["pipeline_path"],
        url=lambda wildcards: config[str(wildcards.SRA)],
    resources:
        mem_mb=4096,
    log:
        "logs/get_nanopore_header_{protocol}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {params.pipeline_path}{output.header})
        samtools view -H {params.url} > {params.pipeline_path}{output.header}
        """


# Body runterladen klappt manuell genau mit diesem Befehl, aber in Snakemake schlaegt es fehl...
rule get_nanopore_body:
    output:
        body="resources/Nanopore/{protocol}/{SRA}/body.sam",
    params:
        pipeline_path=config["pipeline_path"],
        url=lambda wildcards: config[str(wildcards.SRA)],
    resources:
        mem_mb=4096,
    log:
        "logs/get_nanopore_body_{protocol}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view {params.url} | head -n 10000 > {params.pipeline_path}{output.body}
        """


rule combine_nanopore_data:
    input:
        header="resources/Nanopore/{protocol}/{SRA}/header.sam",
        body="resources/Nanopore/{protocol}/{SRA}/body.sam",
    output:
        comb="resources/Nanopore/{protocol}/{SRA}/alignment.sam",
        alignment="resources/Nanopore/{protocol}/{SRA}/alignment.bam",
    params:
        pipeline_path=config["pipeline_path"],
        url=lambda wildcards: config[str(wildcards.SRA)],
    resources:
        mem_mb=4096,
    log:
        "logs/combine_nanopore_data_{protocol}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        cat {params.pipeline_path}{input.header} {params.pipeline_path}{input.body} > {params.pipeline_path}{output.comb}
        samtools view -b {params.pipeline_path}{output.comb} > {params.pipeline_path}{output.alignment}
        """


rule nanopore_bam_index:
    output:
        alignment="resources/Nanopore/{protocol}/{SRA}/alignment.bam.bai",
    params:
        pipeline_path=config["pipeline_path"],
        url=lambda wildcards: config[str(wildcards.SRA)],
    resources:
        mem_mb=4096,
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/nanopore_bam_index_{protocol}_{SRA}.log",
    shell:
        """
        samtools view -b https://ont-open-data.s3.amazonaws.com/gm24385_mod_2021.09/extra_analysis/all.bam | samtools index - {params.pipeline_path}{output}
        """
