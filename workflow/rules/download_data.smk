ref_gene = config.get("sample")
chromosomes = set(chromosome for chromosome in config["seq_platforms"].values())


rule download_genome:
    output:
        "resources/genome.fasta",
    params:
        species=ref_gene["species"],
        datatype=ref_gene["datatype"],
        build=ref_gene["build"],
        release=ref_gene["release"],
    log:
        "logs/data/download_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule genome_index:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/data/genome_index.log",
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
        "logs/data/focus_genome_on_chromosome_{chromosome}.log",
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
        "logs/data/chromosome_index_{chromosome}.log",
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
        "logs/data/rename_chromosome_in_fasta_{chromosome}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_chrom_in_fasta.py"


rule get_fastq_pe:
    output:
        "resources/Illumina_pe/{sample}/{SRA}/{accession}_1.fastq",
        "resources/Illumina_pe/{sample}/{SRA}/{accession}_2.fastq",
    log:
        "logs/data/{sample}/get_fastq_pe_{SRA}_{accession}.log",
    params:
        extra="--skip-technical",
    threads: 6
    # wildcard_constraints:
    #     sample="^(?!simulated_data$).*",
    wrapper:
        "v7.1.0/bio/sra-tools/fasterq-dump"


rule get_fastq_se:
    output:
        "resources/Illumina_se/{sample}/{SRA}/{accession}.fastq",
    log:
        "logs/data/{sample}/get_fastq_se_{SRA}_{accession}.log",
    params:
        extra="--skip-technical",
    threads: 6
    wrapper:
        "v7.1.0/bio/sra-tools/fasterq-dump"


rule trim_fastq_pe:
    input:
        first="resources/Illumina_pe/{sample}/{SRA}/{accession}_1.fastq",
        second="resources/Illumina_pe/{sample}/{SRA}/{accession}_2.fastq",
    output:
        first="resources/Illumina_pe/{sample}/{SRA}/{accession}_1_trimmed.fastq",
        second="resources/Illumina_pe/{sample}/{SRA}/{accession}_2_trimmed.fastq",
    log:
        "logs/data/{sample}/trim_fastq_pe_{SRA}_{accession}.log",
    conda:
        "../envs/fastp.yaml"
    wildcard_constraints:
        sample="^(?!simulated_data$).*",
    shell:
        "fastp --in1 {input.first} --in2 {input.second} --out1 {output.first} --out2 {output.second} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis 2> {log}"


rule trim_fastq_se:
    input:
        first="resources/Illumina_se/{sample}/{SRA}/{accession}.fastq",
    output:
        first="resources/Illumina_se/{sample}/{SRA}/{accession}_trimmed.fastq",
    log:
        "logs/data/{sample}/trim_fastq_se_{SRA}_{accession}.log",
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp --in1 {input.first} --out1 {output.first} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis 2> {log}"


rule get_pacbio_data:
    output:
        alignment="resources/PacBio/{sample}/{SRA}/alignment.bam",
    params:
        url=lambda wildcards: config.get(str(wildcards.SRA)),
        chromosome="chr" + str(config["seq_platforms"].get("PacBio")),
    log:
        "logs/data/{sample}/get_pacbio_data_{SRA}.log",
    resources:
        mem_mb=4096,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -b {params.url} {params.chromosome} > {output.alignment}"
        # "../scripts/get_pacbio_data.py"


# rule get_nanopore_header:
#     output:
#         header="resources/Nanopore/{sample}/{SRA}/header.sam",
#     params:
#
#         url=lambda wildcards: config[str(wildcards.SRA)],
#     resources:
#         mem_mb=4096,
#     log:
#         "logs/get_nanopore_header_{sample}_{SRA}.log",
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         """
#         mkdir -p $(dirname {output.header})
#         samtools view -H {params.url} > {output.header}
#         """


# # Body runterladen klappt manuell genau mit diesem Befehl, aber in Snakemake schlaegt es fehl...
# rule get_nanopore_body:
#     output:
#         body="resources/Nanopore/{sample}/{SRA}/body.sam",
#     params:
#
#         url=lambda wildcards: config[str(wildcards.SRA)],
#     resources:
#         mem_mb=4096,
#     log:
#         "logs/get_nanopore_body_{sample}_{SRA}.log",
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         """
#         samtools view {params.url} | head -n 100000 > {output.body}
#         """


# rule combine_nanopore_data:
#     input:
#         header="resources/Nanopore/{sample}/{SRA}/header.sam",
#         body="resources/Nanopore/{sample}/{SRA}/body.sam",
#     output:
#         comb="resources/Nanopore/{sample}/{SRA}/alignment.sam",
#         alignment="resources/Nanopore/{sample}/{SRA}/alignment.bam",
#     params:
#
#         url=lambda wildcards: config[str(wildcards.SRA)],
#     resources:
#         mem_mb=4096,
#     log:
#         "logs/combine_nanopore_data_{sample}_{SRA}.log",
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         """
#         cat {input.header} {input.body} > {output.comb}
#         samtools view -b {output.comb} > {output.alignment}
#         """


# rule get_nanopore_data:
#     output:
#         "resources/Nanopore/{sample}/{SRA}/alignment.bam",
#     params:
#         url=lambda wildcards: config[str(wildcards.SRA)],
#     resources:
#         mem_mb=4096,
#     log:
#         "logs/data/{sample}/get_nanopore_data_{SRA}.log",
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         """
#         set +o pipefail;
#         mkdir -p $(dirname {output}) 2> {log}
#         samtools view -h {params.url} | head -n 200000 | samtools view -bo {output} - 2> {log}
#         """


# rule nanopore_index:
#     output:
#         alignment="resources/Nanopore/{sample}/{SRA}/alignment.bam.bai",
#     params:
#         url=lambda wildcards: config[str(wildcards.SRA)],
#     resources:
#         mem_mb=4096,
#     conda:
#         "../envs/samtools.yaml"
#     log:
#         "logs/data/{sample}/nanopore_index_{SRA}.log",
#     shell:
#         "samtools view -b {params.url} | samtools index - {output} 2> {log}"


# TODO: Does not work for replicate2. You have to download this manually with wget right now
rule get_nanopore_data:
    output:
        alignment="resources/Nanopore/{sample}/{SRA}/alignment.bam",
    params:
        url=lambda wildcards: config.get(str(wildcards.SRA)),
        chromosome="chr" + str(config["seq_platforms"].get("Nanopore")),
    log:
        "logs/data/{sample}/get_nanopore_data_{SRA}.log",
    resources:
        mem_mb=4096,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -b {params.url} {params.chromosome} > {output.alignment}"
