chr_chromosome = "chr" + chromosome_conf["chromosome"]


rule download_bedGraphs:
    output:
        "resources/HG002/{bedGraph}.bedGraph.gz",
    log:
        "logs/download_bedGraphs{bedGraph}.log",
    params:
        pipeline_path=config["pipeline_path"],
        bedGraphs = config["bedGraphs_HG002"]
    script:
        "../scripts/get_bedGraph_data.py"

ruleorder: filter_bedGraphs > extract_data

rule extract_data:
    input:
        "resources/HG002/{bedGraph}.bedGraph.gz"
    output:
        "resources/HG002/{bedGraph}.bedGraph"
    log:
        "logs/process_data{bedGraph}.log",
    shell:
        "gunzip -c {input} > {output}"

rule filter_bedGraphs:
    input:
        "resources/HG002/{bedGraph}.bedGraph"
    output:
        "resources/HG002/{bedGraph}_filtered.bedGraph"
    log:
        "logs/filter_bedGraphs{bedGraph}.log",
    params:
        chromosome=chr_chromosome
    shell:
        """
        awk '$1 == "{params.chromosome}" {{print}}' {input} > {output}
        """

rule compute_avg_bedGraph:
    input:
        # expand("resources/{SRA}/HG002/bedGraph/{bed}.bedGraph", bed=bedGraphs),
        expand("resources/HG002/{bedGraph}_filtered.bedGraph", bedGraph=config["bedGraphs_HG002"])

    output:
        "resources/bed_avg.bedGraph",
    log:
        "logs/compute_avg_bedGraph.log",
    params:
        chromosome=chr_chromosome,
    script:
        "../scripts/compute_avg_bedGraph.py"

rule plot_avg_bedGraph:
    input:
        candidates="resources/candidates.vcf",
        bedgraphs=expand("resources/HG002/{bedGraph}_filtered.bedGraph", bedGraph=config["bedGraphs_HG002"])
    output:
        cov="resources/bed_avg_cov.png",
        meth="resources/bed_avg_meth.png",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_avg_bedGraph.log",
    script:
        "../scripts/plot_avg_bedGraph.py"
