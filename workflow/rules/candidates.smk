rule find_candidates:
    input:
        varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
        fasta="resources/chromosome_{chromosome}.fasta",
    output:
        "resources/{chromosome}/candidates.bcf",
    log:
        "logs/candidates/find_candidates_{chromosome}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} methylation-candidates {input.fasta} {output} 2> {log}"


rule split_candidates:
    input:
        "resources/{chrom}/candidates.bcf",
    output:
        scatter.split_candidates("resources/{{chrom}}/candidates_{scatteritem}.bcf"),
    log:
        "logs/candidates/split_candidates_{chrom}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output} 2> {log}"
