# rule compute_meth_observations:
#     input:
#         chromosome=lambda wildcards: expand(
#             "resources/chromosome_{chrom}.fasta",
#             chrom=chromosome_by_platform[wildcards.platform],
#         ),
#         genome_index=lambda wildcards: expand(
#             "resources/chromosome_{chrom}.fasta.fai",
#             chrom=chromosome_by_platform[wildcards.platform],
#         ),
#         alignments="resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam",
#         alignment_index="resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam.bai",
#         candidates=lambda wildcards: expand(
#             "resources/{{platform}}/{{protocol}}/21/candidates_{{scatteritem}}.bcf",
#             chrom=chromosome_by_platform[wildcards.platform],
#         ),
#     output:
#         "results/{platform}/{protocol}/normal_{scatteritem}.bcf",
#     log:
#         "logs/compute_meth_observations_{platform}_{protocol}_{scatteritem}.log",
#     conda:
#         "../envs/varlociraptor.yaml"
#     params:
#         varlo_path=config["varlo_path"],
#         pipeline_path=config["pipeline_path"],
#     shell:
#         """
#         cd {params.varlo_path}
#         if [[ "{wildcards.platform}" == "Illumina_pe" || "{wildcards.platform}" == "Illumina_se" ]]; then
#             PLATFORM="Illumina"
#         else
#             PLATFORM="{wildcards.platform}"
#         fi
#         echo $PLATFORM
#         cargo run --release -- preprocess variants {params.pipeline_path}{input.chromosome} --candidates {params.pipeline_path}{input.candidates} --bam {params.pipeline_path}{input.alignments} --read-type $PLATFORM --max-depth 5000 > {params.pipeline_path}{output}
#         """


rule compute_meth_observations:
    input:
        chromosome=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
        genome_index=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
        alignments="resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam",
        alignment_index="resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam.bai",
#         candidates=lambda wildcards: expand(
#             "resources/{{platform}}/{{protocol}}/21/candidates_{{scatteritem}}.bcf",
#             chrom=chromosome_by_platform[wildcards.platform],
#         ),
        candidates=lambda wildcards: expand(
            "resources/{{platform}}/{{protocol}}/{chrom}/candidates_{{scatteritem}}.bcf",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
    output:
        "results/{platform}/{protocol}/normal_{scatteritem}.bcf",
    log:
        "logs/compute_meth_observations_{platform}_{protocol}_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        varlo_path=config["varlo_path"],
        pipeline_path=config["pipeline_path"],
    shell:
        """
        cd {params.varlo_path}
        if [[ "{wildcards.platform}" == "Illumina_pe" || "{wildcards.platform}" == "Illumina_se" ]]; then
            PLATFORM="Illumina"
        else
            PLATFORM="{wildcards.platform}"
        fi
        echo $PLATFORM
        cargo run --release -- preprocess variants {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --read-type $PLATFORM --max-depth 5000 > {output}
        """
        # cargo run --release -- preprocess variants --omit-mapq-adjustment {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --read-type $PLATFORM --max-depth 5000 > {output}


rule call_methylation:
    input:
        preprocess_obs="results/{platform}/{protocol}/normal_{scatteritem}.bcf",
        scenario="resources/scenario.yaml",
    output:
        "results/{platform}/{protocol}/calls_{scatteritem}.bcf",
    log:
        "logs/call_methylation_{platform}_{protocol}_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        varlo_path=config["varlo_path"],
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        cd {params.varlo_path}
        cargo run --release -- call variants generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output}
        """


# TODO: Reactivate, right now it deletes too much data
rule filter_calls:
    input:
        "results/{platform}/{protocol}/calls_{scatteritem}.bcf",
    output:
        "results/{platform}/{protocol}/calls_{scatteritem}.filtered.bcf",
    log:
        "logs/filter_calls_{platform}_{protocol}_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        varlo_path=config["varlo_path"],
        pipeline_path=config["pipeline_path"],
        event="PRESENT",
    shell:
        """
        cd {params.varlo_path}
        cargo run --release -- filter-calls control-fdr --mode local-smart {input} --events {params.event} --fdr 0.005 > {output}
        """


rule calls_to_vcf:
    input:
        "results/{platform}/{protocol}/calls_{scatteritem}.bcf",
    output:
        "results/{platform}/{protocol}/calls_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/convert_to_vcf_{platform}_{protocol}_{scatteritem}.log",
    threads: 10
    shell:
        """
        bcftools view --threads {threads} {input} -o {output}
        """


rule gather_calls:
    input:
        gather.split_candidates(
            "results/{{platform}}/{{protocol}}/calls_{scatteritem}.vcf"
        ),
    output:
        "results/{platform}/{protocol}/varlo.vcf",
    log:
        "logs/gather_calls_{platform}_{protocol}.log",
    conda:
        "../envs/cat.yaml"
    shell:
        "cat {input} > {output}"


rule rename_varlo_output:
    input:
        "results/{platform}/{protocol}/varlo.vcf",
    output:
        "results/{platform}/{protocol}/result_files/varlo.bed",
    shell:
        "mv {input} {output}"
