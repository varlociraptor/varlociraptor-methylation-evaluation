# Use Varlociraptor to find methylation candidates in the reference genome
rule find_candidates:
    input:
        "resources/chromosome_{chromosome}.fasta",
    output:
        "resources/{chromosome}/candidates.bcf",
    log:
        "logs/candidates/find_candidates/{chromosome}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor methylation-candidates {input} --motif CG {output} 2> {log}"


rule split_candidates:
    input:
        "resources/{chrom}/candidates.bcf",
    output:
        scatter.split_candidates("resources/{{chrom}}/candidates_{scatteritem}.bcf"),
    log:
        "logs/candidates/split_candidates/{chrom}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output} 2> {log}"


rule index_candidates:
    input:
        "resources/{chrom}/candidates.bcf",
    output:
        "resources/{chrom}/candidates.bcf.csi",
    log:
        "logs/candidates/index_candidates/{chrom}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools index {input}"
