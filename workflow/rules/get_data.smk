chromosome_conf = config["sample"]

rule get_genome:
    output:
        "resources/genome.fasta",
    params:
        species=chromosome_conf["species"],
        datatype=chromosome_conf["datatype"],
        build=chromosome_conf["build"],
        release=chromosome_conf["release"],
    log:
        "logs/get_chromosome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule get_chromosome:
    output:
        "resources/chrom.fasta",
    params:
        species=chromosome_conf["species"],
        datatype=chromosome_conf["datatype"],
        build=chromosome_conf["build"],
        release=chromosome_conf["release"],
        chromosome=chromosome_conf["chromosome"],
    log:
        "logs/get_chromosome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule filter_genome:
    input:
        "resources/genome.fasta",
    output:
        "resources/chromosome.fasta",
    log:
        "logs/filter_genome.log",
    conda:
        "../envs/samtools.yaml"
    params:
        chromosome=chromosome_conf["chromosome"],
    threads: 10
    shell:
        """ 
        samtools faidx {input} {params.chromosome} > {output}
        """

rule genome_index:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/chromosome_index.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        samtools faidx {params.pipeline_path}{input}
        """

rule get_pacbio_data:
    output:
        alignment="resources/PacBio/{protocol}/{SRA}/alignment.bam"
    params:
        url=lambda wildcards: config[str(wildcards.SRA)]
    resources: mem_mb=4096
    script:
        "../scripts/get_pacbio_data.py"

rule nanopore_bam_index:
    output:
        alignment="resources/Nanopore/{protocol}/{SRA}/alignment.bam.bai"
    params:
        pipeline_path=config["pipeline_path"],
        url=lambda wildcards: config[str(wildcards.SRA)]
    resources: mem_mb=4096
    conda: 
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -b https://ont-open-data.s3.amazonaws.com/gm24385_mod_2021.09/extra_analysis/all.bam | samtools index - {params.pipeline_path}{output}
	    """

# Body runterladen klappt manuell genau mit diesem Befehl, aber in Snakemake schlaegt es fehl...
rule get_nanopore_data:
    output:
        header="resources/Nanopore/{protocol}/{SRA}/header.sam",
        body="resources/Nanopore/{protocol}/{SRA}/body.sam",
    params:
        pipeline_path=config["pipeline_path"],
        url=lambda wildcards: config[str(wildcards.SRA)]
    resources: mem_mb=4096
    conda: 
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {params.pipeline_path}{output.header})
        samtools view -H {params.url} > {params.pipeline_path}{output.header}
        samtools view {params.url} | head -n 1000 > {params.pipeline_path}{output.body}
	    """

rule combine_nanopore_data:
    input:
        header="resources/Nanopore/{protocol}/{SRA}/header.sam",
        body="resources/Nanopore/{protocol}/{SRA}/body.sam",
    output:
        comb="resources/Nanopore/{protocol}/{SRA}/combined.sam",
        alignment="resources/Nanopore/{protocol}/{SRA}/alignment.bam"
    params:
        pipeline_path=config["pipeline_path"],
        url=lambda wildcards: config[str(wildcards.SRA)]
    resources: mem_mb=4096
    conda: 
        "../envs/samtools.yaml"
    shell:
        """
        cat {params.pipeline_path}{input.header} {params.pipeline_path}{input.body} > {params.pipeline_path}{output.comb}
        samtools view -b {params.pipeline_path}{output.comb} > {params.pipeline_path}{output.alignment}
	    """
        # samtools view {params.url} | head -n 10000 | samtools view -b > {params.pipeline_path}{output.alignment}
        # samtools view -h {params.url} chr21 -O bam > {params.pipeline_path}{output}


rule get_fastq_pe:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "resources/Illumina/{protocol}/{SRA}/{accession}_1.fastq",
        "resources/Illumina/{protocol}/{SRA}/{accession}_2.fastq",
    log:
        "logs/pe/{accession}{protocol}{SRA}.log"
    params:
        extra="--skip-technical"
    threads: 6  # defaults to 6
    # conda:
    #     "../envs/fastq-wrapper.yaml"
    wrapper:
        "v3.0.2/bio/sra-tools/fasterq-dump"



rule trim_fastq_pe:
    input:
        first="resources/Illumina/{protocol}/{SRA}/{accession}_1.fastq",
        second="resources/Illumina/{protocol}/{SRA}/{accession}_2.fastq",
    output:
        first="resources/Illumina/{protocol}/{SRA}/{accession}_1_trimmed.fastq",
        second="resources/Illumina/{protocol}/{SRA}/{accession}_2_trimmed.fastq",
    log:
        "logs/trim_fastq_pe_{protocol}{SRA}_{accession}.log",
    conda:
        "../envs/fastp.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        fastp --in1 {input.first} --in2 {input.second} --out1 {output.first} --out2 {output.second} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis
        """


  
