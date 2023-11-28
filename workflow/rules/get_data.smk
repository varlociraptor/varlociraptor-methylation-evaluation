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
        "resources/pacbio/HG002.GRCh38.haplotagged.bam",
    log:
        "logs/get_pacbio_data.log",
    script:
        "../scripts/get_pacbio_data.py"

rule get_fastq_pe:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "resources/{platform}/{SRA}/{accession}_1.fastq",
        "resources/{platform}/{SRA}/{accession}_2.fastq",
    log:
        "logs/pe/{accession}{platform}{SRA}.log"
    params:
        extra="--skip-technical"
    threads: 6  # defaults to 6
    conda:
        "../envs/fastq-wrapper.yaml"
    wrapper:
        "v2.6.0/bio/sra-tools/fasterq-dump"


rule trim_fastq_pe:
    input:
        first="resources/{platform}/{SRA}/{accession}_1.fastq",
        second="resources/{platform}/{SRA}/{accession}_2.fastq",
    output:
        first="resources/{platform}/{SRA}/{accession}_1_trimmed.fastq",
        second="resources/{platform}/{SRA}/{accession}_2_trimmed.fastq",
    log:
        "logs/trim_fastq_pe_{platform}{SRA}_{accession}.log",
    conda:
        "../envs/fastp.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        fastp --in1 {input.first} --in2 {input.second} --out1 {output.first} --out2 {output.second} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis
        """


  