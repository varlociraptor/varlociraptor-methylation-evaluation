rule download_bedGraph:
    output:
        "resources/HG002/{bedGraph}.bedGraph.gz",
    log:
        "logs/avg_bedgraph/download_bedGraph_{bedGraph}.log",
    conda:
        "../envs/python.yaml"
    params:
        bedGraphs=config["bedGraphs_HG002"],
        bedgraph_path="resources/HG002",
    script:
        "../scripts/get_bedGraph_data.py"


ruleorder: filter_bedGraph_on_chromosome > unzip_bedgraph


rule unzip_bedgraph:
    input:
        "resources/HG002/{bedGraph}.bedGraph.gz",
    output:
        "resources/HG002/{bedGraph}.bedGraph",
    log:
        "logs/avg_bedgraph/unzip_bedgraph_{bedGraph}.log",
    conda:
        "../envs/shell_cmds.yaml"
    shell:
        "gunzip -c {input} > {output} 2> {log}"


rule filter_bedGraph_on_chromosome:
    input:
        "resources/HG002/{bedGraph}.bedGraph",
    output:
        "resources/HG002/{bedGraph}-{chromosome}.bedGraph",
    wildcard_constraints:
        chromosome="[^_]+",
    log:
        "logs/avg_bedgraph/filter_bedGraph_{bedGraph}_on_chromosome_{chromosome}.log",
    conda:
        "../envs/shell_cmds.yaml"
    shell:
        "awk '$1 == \"{wildcards.chromosome}\" {{print}}' {input} > {output} 2> {log}"


# Strategy: Go through all 16 bedgraphs.
# We compute the average methylation level for a position if this position appears in at least 12 of the bedgraphs and only if the maximum difference in methylation levels is <= 20 between all begraphs.
# We compute the average methylation level by dividing all methylated reads by the total number of reads in all bedgraphs at this position
rule compute_avg_bedGraph:
    input:
        expand(
            "resources/HG002/{bedGraph}-{chromosome}.bedGraph",
            bedGraph=config["bedGraphs_HG002"],
            chromosome=[chrom for chrom in chr_chromosome],
        ),
    output:
        "resources/bed_avg_{chromosome}.bedGraph",
    wildcard_constraints:
        chromosome="[^_]+",
    log:
        "logs/avg_bedgraph/compute_avg_bedGraph_{chromosome}.log",
    conda:
        "../envs/python.yaml"
    params:
        chromosome=lambda wildcards: wildcards.chromosome,
    script:
        "../scripts/compute_avg_bedGraph.py"
