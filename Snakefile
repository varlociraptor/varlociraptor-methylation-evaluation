rule get_genome:
    output:
        "resources/genome.fasta",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="110",
    log:
        "logs/get_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"

rule get_chromosome:
    output:
        "resources/chr1.fasta",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="110",
        chromosome="I",  # optional: restrict to chromosome
        # branch="plants",  # optional: specify branch
    log:
        "logs/get_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"

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