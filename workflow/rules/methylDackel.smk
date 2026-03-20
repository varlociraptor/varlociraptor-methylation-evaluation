# It does not work with chromosome-wise fasta files, so we use the genome
rule methylDackel_compute_meth:
    input:
        genome=lambda wildcards: (
            expand(
                "resources/chromosome_{chrom}.fasta",
                chrom=config["seq_platforms"].get(wildcards.platform),
            )
            if wildcards.sample.startswith("simulated_data")
            else ["resources/genome.fasta"]
        ),
        genome_index=lambda wildcards: (
            expand(
                "resources/chromosome_{chrom}.fasta.fai",
                chrom=config["seq_platforms"].get(wildcards.platform),
            )
            if wildcards.sample.startswith("simulated_data")
            else ["resources/genome.fasta.fai"]
        ),
        alignment="resources/{platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam",
        alignment_index="resources/{platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam.bai",
    output:
        temp("results/single_sample/{platform}/called/{sample}/result_files/alignments_CpG_{scatteritem}.bedGraph"),
    conda:
        "../envs/methylDackel.yaml"
    log:
        "logs/methylDackel/methylDackel_compute_meth/{platform}_{sample}_{scatteritem}.log",
    benchmark:
        repeat("benchmarks/{platform}/methylDackel/methylDackel_compute_meth/{sample}_{scatteritem}.bwa.benchmark.txt", 3)
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output})
        OUTDIR=$(dirname {output})/alignments_{wildcards.scatteritem}
        mkdir -p "$OUTDIR"
        MethylDackel extract {input.genome} {input.alignment} -o "$OUTDIR" --mergeContext 2> {log}
        mv "$OUTDIR"_CpG.bedGraph {output}
        """


rule methylDackel_gather_meth:
    input:
        gather.split_candidates(
            "results/single_sample/{{platform}}/called/{{sample}}/result_files/alignments_CpG_{scatteritem}.bedGraph"
        ),
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/alignments_CpG.bedGraph",
    log:
        "logs/methylDackel/methylDackel_gather_meth/{platform}_{sample}.log",
    conda:
        "../envs/general.yaml"
    benchmark:
        repeat("benchmarks/{platform}/methylDackel/methylDackel_gather_meth/{platform}_{sample}.log", 3)
    shell:
        """
        head -n1 $(echo {input} | tr ' ' '\n' | head -n1) > {output} 2> {log}
        for f in {input}; do tail -n +2 "$f"; done >> {output} 2>> {log}
        """


rule methylDackel_rename_output:
    input:
        "results/single_sample/{platform}/called/{sample}/result_files/alignments_CpG.bedGraph",
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/methylDackel.bed",
    log:
        "logs/methylDackel/methylDackel_rename_output/{platform}_{sample}.log",
    conda:
        "../envs/general.yaml"
    shell:
        "mkdir -p $(dirname {log}) && mkdir -p $(dirname {output}) && mv {input} {output} 2> {log}"
