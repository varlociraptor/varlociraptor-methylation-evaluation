configfile: "config/config.yaml"


sample_config = config["sample"]
resource_path = config["resource_path"]


rule get_genome:
    output:
        "resources/genome.fasta",
    params:
        species=sample_config["species"],
        datatype=sample_config["datatype"],
        build=sample_config["build"],
        release=sample_config["release"],
    log:
        "logs/get_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule get_chromosome:
    output:
        "resources/chr1.fasta",
    params:
        species=sample_config["species"],
        datatype=sample_config["datatype"],
        build=sample_config["build"],
        release=sample_config["release"],
        chromosome=sample_config["chromosome"],  # optional: restrict to chromosome
        # branch="plants",  # optional: specify branch
    log:
        "logs/get_chromosome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule find_candidates:
    input:
        "resources/example-genome.fasta",
    output:
        "resources/example-candidates.bcf",
    log:
        "logs/find_candidates.log",
    shell:
        """ 
        cd ~/Documents/Promotion/varlociraptor/
        cargo run -- methylation-candidates {resource_path}/{input} {resource_path}/{output}
        """


rule sam_to_bam:
    input:
        "resources/example-reads.sam",
    output:
        "resources/example-reads.bam",
    log:
        "logs/sam_to_bam.log",
    shell:
        "samtools view -bS resources/example-reads.sam > resources/example-reads.bam"


rule compute_meth_probs:
    input:
        fasta="resources/example-genome.fasta",
        bam="resources/example-reads.bam",
        candidates="resources/example-candidates.bcf",
    output:
        "resources/observations.bcf",
    log:
        "logs/copute.log",
    shell:
        """ 
        cd ~/Documents/Promotion/varlociraptor/
        cargo run -- preprocess variants {resource_path}/{input.fasta} --candidates {resource_path}/{input.candidates} --bam {resource_path}/{input.bam} > {resource_path}/{output}
        """
