configfile: "config/config.yaml"

rule get_genome:
    output:
        "resources/genome.fasta",
    params:
        species=config["sample"]["species"],
        datatype=config["sample"]["datatype"],
        build=config["sample"]["build"],
        release=config["sample"]["release"],
    log:
        "logs/get_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"

rule get_chromosome:
    output:
        "resources/chr1.fasta",
    params:
        species=config["sample"]["species"],
        datatype=config["sample"]["datatype"],
        build=config["sample"]["build"],
        release=config["sample"]["release"],
        chromosome=config["sample"]["chromosome"],  # optional: restrict to chromosome
        # branch="plants",  # optional: specify branch
    log:
        "logs/get_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


"""Later    shell: varlociraptor --candidates {input} {output}"""

rule find_methylation:
    input: 
        "resources/genome.fasta",
    output:
        "resources/candidates.vcf",
    log:
        "logs/compute_candidates.log",
    shell: 
        """ 
        cd ~/Documents/Promotion/varlociraptor/
        cargo run -- candidates ~/Documents/Promotion/varlociraptor-methylation-evaluation/{input} ~/Documents/Promotion/varlociraptor-methylation-evaluation/{output}
        """
        

