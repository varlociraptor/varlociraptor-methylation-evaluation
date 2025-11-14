rule varlociraptor_preprocess:
    input:
        chromosome=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=chromosome_by_seq_platform.get(wildcards.seq_platform),
        ),
        genome_index=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=chromosome_by_seq_platform.get(wildcards.seq_platform),
        ),
        alignments="resources/{seq_platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam",
        alignment_index="resources/{seq_platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam.bai",
        candidates=lambda wildcards: expand(
            "resources/{chrom}/candidates_{{scatteritem}}.bcf",
            chrom=chromosome_by_seq_platform.get(wildcards.seq_platform),
        ),
    output:
        "results/preprocessed/{seq_platform}/{sample}/normal_{scatteritem}.bcf",
    log:
        "logs/varlociraptor_single/varlociraptor_preprocess/{seq_platform}_{sample}_{scatteritem}.log",
    benchmark:
        "benchmarks/{seq_platform}/varlociraptor/preprocessing/{sample}_{scatteritem}.bwa.benchmark.txt"
    conda:
        "../envs/varlociraptor.yaml"
    resources:
        mem_mb=16000,
    shell:
        """
        if [[ "{wildcards.seq_platform}" == "Illumina_pe" || "{wildcards.seq_platform}" == "Illumina_se" ]]; then
            meth="converted"
        else
            meth="annotated"
        fi
        varlociraptor preprocess variants --omit-mapq-adjustment {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --methylation-read-type $meth --max-depth 5000 > {output} 2> {log}
        """


rule varlociraptor_call:
    input:
        preprocess_obs="results/preprocessed/{seq_platform}/{sample}/normal_{scatteritem}.bcf",
        scenario="resources/scenarios/scenario.yaml",
    output:
        "results/single_sample/{seq_platform}/called/{sample}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor_single/varlociraptor_call/{seq_platform}_{sample}_{scatteritem}.log",
    benchmark:
        "benchmarks/{seq_platform}/varlociraptor/calling/{sample}_{scatteritem}.bwa.benchmark.txt"
    conda:
        "../envs/varlociraptor.yaml"
    wildcard_constraints:
        seq_platform="(?!multi_sample).*",
    shell:
        "varlociraptor call variants generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output} 2> {log}"


# TODO: Skip this step, right now it would be useless since I debug so much
rule calls_to_vcf:
    input:
        "results/{call_type}/{seq_platform}/{sample}/calls_{scatteritem}.bcf",
    output:
        "results/{call_type}/{seq_platform}/{sample}/calls_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/varlociraptor_single/calls_to_vcf/{call_type}_{seq_platform}_{sample}_{scatteritem}.log",
    threads: 10
    shell:
        # "touch {output} 2> {log}"
        "bcftools view --threads {threads} {input} -o {output} 2> {log}"


rule gather_calls:
    input:
        gather.split_candidates(
            "results/{{call_type}}/{{seq_platform}}/{{sample}}/calls_{scatteritem}.vcf"
        ),
    output:
        "results/{call_type}/{seq_platform}/{sample}/result_files/varlo.bed",
    log:
        "logs/varlociraptor_single/gather_calls/{call_type}_{seq_platform}_{sample}.log",
    conda:
        "../envs/general.yaml"
    shell:
        "cat {input} > {output} 2> {log}"


# rule rename_varlo_output:
#     input:
#         "results/{seq_platform}/{sample}/varlo.vcf",
#     output:
#         "results/{seq_platform}/{sample}/result_files/varlo.bed",
#     log:
#         "logs/varlociraptor/{seq_platform}/{sample}/rename_varlo_output.log",
#     shell:
#         "mv {input} {output}"
