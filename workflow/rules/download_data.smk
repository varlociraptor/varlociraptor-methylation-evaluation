ref_gene = config.get("sample", {})
chromosomes = set(chromosome for chromosome in config["seq_platforms"].values())


rule bcf_to_vcf:
    input:
        "{file}.bcf",
    output:
        "{file}.vcf",
    log:
        "logs/bcf_to_vcf/{file}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        # "touch {output} 2> {log}"
        "bcftools view --threads {threads} {input} -o {output} 2> {log}"


rule download_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/download_data/download_genome/download.log",
    cache: "omit-software"
    params:
        species=ref_gene.get("species"),
        datatype=ref_gene.get("datatype"),
        build=ref_gene.get("build"),
        release=ref_gene.get("release"),
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
        "resources/{chromosome}.fasta",
    log:
        "logs/download_data/focus_genome_on_chromosome/{chromosome}.log",
    wildcard_constraints:
        chromosome="(?!.*genome$)[^/]+",
    conda:
        "../envs/samtools.yaml"
    threads: 6
    shell:
        """
        if [[ {wildcards.chromosome} == genome ]]; then \
        echo "Copying genome fasta to output"
        cp {input} {output}
        else
            echo "Extracting chromosome {wildcards.chromosome} from genome fasta"
            samtools faidx {input} {wildcards.chromosome} > {output}
        fi 2> {log}
        """


rule unzip_genome:
    input:
        "resources/{genome}.fasta",
    output:
        "resources/{genome}.fasta",
    log:
        "logs/download_data/unzip_genome/{genome}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "gunzip -c {input} > {output} 2> {log}"


rule chromosome_index:
    input:
        "resources/{chromosome}.fasta",
    output:
        "resources/{chromosome}.fasta.fai",
    log:
        "logs/download_data/chromosome_index/{chromosome}.log",
    wildcard_constraints:
        chromosome="(?!.*genome$)[^/]+",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input} 2> {log}"


rule rename_chromosome_in_fasta:
    input:
        "resources/{chromosome}.fasta",
    output:
        "resources/chr_{chromosome}.fasta",
    log:
        "logs/download_data/rename_chromosome_in_fasta/{chromosome}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_chrom_in_fasta.py"


# We need to call the wildcard accession and not SRA because of the wrapper
rule get_fastq_pe:
    output:
        "resources/Illumina_pe/{sample}/{SRA}/{accession}_1.fastq.gz",
        "resources/Illumina_pe/{sample}/{SRA}/{accession}_2.fastq.gz",
    log:
        "logs/download_data/get_fastq_pe/{sample}_{SRA}_{accession}.log",
    threads: 6
    params:
        extra="--skip-technical",
    # wildcard_constraints:
    #     sample="^(?!simulated_data).*",
    wrapper:
        "v7.1.0/bio/sra-tools/fasterq-dump"


rule get_fastq_se:
    output:
        "resources/Illumina_se/{sample}/{SRA}/{accession}.fastq.gz",
    log:
        "logs/download_data/get_fastq_se/{sample}_{SRA}_{accession}.log",
    threads: 6
    params:
        extra="--skip-technical",
    wrapper:
        "v7.1.0/bio/sra-tools/fasterq-dump"


rule trim_fastq_pe:
    input:
        sample=[
            "resources/Illumina_pe/{sample}/{SRA}/{SRA}_1.fastq.gz",
            "resources/Illumina_pe/{sample}/{SRA}/{SRA}_2.fastq.gz",
        ],
    output:
        trimmed=[
            "resources/Illumina_pe/{sample}/{SRA}/{SRA}_1_trimmed.fastq.gz",
            "resources/Illumina_pe/{sample}/{SRA}/{SRA}_2_trimmed.fastq.gz",
        ],
        merged="trimmed/pe/{sample}_{SRA}.merged.fastq.gz",
        failed="trimmed/pe/{sample}_{SRA}.failed.fastq.gz",
        html="report/pe/{sample}_{SRA}.html",
        json="report/pe/{sample}_{SRA}.json",
    log:
        "logs/fastp/pe/{sample}_{SRA}.log",
    threads: 8
    params:
        # adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge",
    wrapper:
        "v9.4.1/bio/fastp"


rule trim_fastq_se:
    input:
        first="resources/Illumina_se/{sample}/{SRA}/{accession}.fastq.gz",
    output:
        first="resources/Illumina_se/{sample}/{SRA}/{accession}_trimmed.fastq.gz",
    log:
        "logs/download_data/trim_fastq_se/{sample}_{SRA}_{accession}.log",
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp --in1 {input.first} --out1 {output.first} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis 2> {log}"


rule get_pacbio_data:
    output:
        alignment="resources/PacBio/{sample}/{SRA}/alignment.bam",
    log:
        "logs/download_data/get_pacbio_data/{sample}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=4096,
    params:
        url=lambda wildcards: config.get(str(wildcards.SRA)),
        chromosome=f"chr{config['seq_platforms'].get('PacBio')}",
    shell:
        "samtools view -b {params.url} {params.chromosome} > {output.alignment} 2> {log}"


rule get_nanopore_index:
    output:
        "resources/Nanopore/{sample}/{SRA}/alignment.bam.bai",
    log:
        "logs/download_data/get_nanopore_index/{sample}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        url=lambda wildcards: config.get(str(wildcards.SRA)),
    shell:
        "wget {params.url}.bai -O {output} 2> {log}"


# TODO: Does not work for replicate2. You have to download this manually with wget right now
rule get_nanopore_data:
    output:
        alignment="resources/Nanopore/{sample}/{SRA}/alignment.bam",
    log:
        "logs/download_data/get_nanopore_data/{sample}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=4096,
    params:
        url=lambda wc: config.get(str(wc.SRA)),
        chromosome=lambda wc: f"chr{config['seq_platforms']['Nanopore']}",
    shell:
        """
        mkdir -p $(dirname {output.alignment}) \
         && wget -qO- {params.url} \
        | samtools view -b - {params.chromosome} > {output.alignment} 2> {log}
        """
