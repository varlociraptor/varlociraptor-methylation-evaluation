chromosome_conf = config["sample"]

rule get_genome:
    output:
        "resources/{SRA}/genome.fasta",
    params:
        species=chromosome_conf["species"],
        datatype=chromosome_conf["datatype"],
        build=chromosome_conf["build"],
        release=chromosome_conf["release"],
    log:
        "logs/get_chromosome{SRA}.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule get_chromosome:
    output:
        "resources/{SRA}/chrom.fasta",
    params:
        species=chromosome_conf["species"],
        datatype=chromosome_conf["datatype"],
        build=chromosome_conf["build"],
        release=chromosome_conf["release"],
        chromosome=chromosome_conf["chromosome"],
    log:
        "logs/get_chromosome{SRA}.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule filter_genome:
    input:
        "resources/{SRA}/genome.fasta",
    output:
        "resources/{SRA}/chromosome.fasta",
    log:
        "logs/filter_genome{SRA}.log",
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
        "resources/{SRA}/genome.fasta",
    output:
        "resources/{SRA}/genome.fasta.fai",
    log:
        "logs/chromosome_index{SRA}.log",
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
        # the wildcard name must be accession, pointing to an SRA number
        "resources/{SRA}/{accession}_1.fastq",
        "resources/{SRA}/{accession}_2.fastq",
    log:
        "logs/pe/{accession}{SRA}.log"
    params:
        extra="--skip-technical"
    threads: 6  # defaults to 6
    conda:
        "../envs/fastq-wrapper.yaml"
    wrapper:
        "v2.6.0/bio/sra-tools/fasterq-dump"