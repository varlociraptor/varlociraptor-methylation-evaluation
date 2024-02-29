chr_chromosome = "chr" + chromosome_conf["chromosome"]


rule download_bedGraphs:
    output:
        "resources/HG002/{bedGraph}.bedGraph.gz",
    log:
        "logs/download_bedGraphs_{bedGraph}.log",
    conda:
        "../envs/python.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        bedGraphs=config["bedGraphs_HG002"],
        bedgraph_path="resources/HG002",
    script:
        "../scripts/get_bedGraph_data.py"


ruleorder: filter_bedGraphs > extract_data


rule extract_data:
    input:
        "resources/HG002/{bedGraph}.bedGraph.gz",
    output:
        "resources/HG002/{bedGraph}.bedGraph",
    log:
        "logs/extract_data_{bedGraph}.log",
    conda:
        "../envs/gunzip.yaml"
    shell:
        "gunzip -c {input} > {output}"


rule filter_bedGraphs:
    input:
        "resources/HG002/{bedGraph}.bedGraph",
    output:
        "resources/HG002/{bedGraph}-{chromosome}.bedGraph",
    wildcard_constraints:
        chromosome="[^_]+",
    log:
        "logs/filter_bedGraphs_{bedGraph}_{chromosome}.log",
    conda:
        "../envs/awk.yaml"
    params:
        chromosome=chr_chromosome,
    shell:
        """
        awk '$1 == "{params.chromosome}" {{print}}' {input} > {output}
        """


rule compute_avg_bedGraph:
    input:
        # expand("resources/{SRA}/HG002/bedGraph/{bed}.bedGraph", bed=bedGraphs),
        expand(
            "resources/HG002/{bedGraph}-{chromosome}.bedGraph",
            bedGraph=config["bedGraphs_HG002"],
            chromosome=chr_chromosome,
        ),
    output:
        "resources/bed_avg_{chromosome}.bedGraph",
    wildcard_constraints:
        chromosome="[^_]+",
    log:
        "logs/compute_avg_bedGraph_{chromosome}.log",
    conda:
        "../envs/python.yaml"
    params:
        chromosome=chr_chromosome,
    script:
        "../scripts/compute_avg_bedGraph.py"


rule plot_avg_bedGraph:
    input:
        candidates=expand(
            "resources/{chro}/candidates.vcf", chro=chromosome_conf["chromosome"]
        ),
        bedgraphs=expand(
            "resources/HG002/{bedGraph}-{chromosome}.bedGraph",
            bedGraph=config["bedGraphs_HG002"],
            chromosome=chr_chromosome,
        ),
    output:
        cov="resources/bed_avg_cov.png",
        meth="resources/bed_avg_meth.png",
    conda:
        "../envs/plot.yaml"
    wildcard_constraints:
        chromosome="[^_]+",
    log:
        "logs/plot_avg_bedGraph.log",
    script:
        "../scripts/plot_avg_bedGraph.py"
