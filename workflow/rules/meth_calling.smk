# TODO: Missingoutputexception
rule download_varlociraptor:
    output:
        directory("resources/tools/varlociraptor"),
    log:
        "../logs/download_varlociraptor.log",
    conda:
        "../envs/sheel_cmds.yaml"
    shell:
        """
        mkdir -p resources/tools
        cd resources/tools
        git clone https://github.com/varlociraptor/varlociraptor.git
        cd varlociraptor
        git checkout methylation-paired-end-master-new
        """


rule compute_meth_observations:
    input:
        varlo=directory("resources/tools/varlociraptor"),
        chromosome=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
        genome_index=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
        alignments="resources/{platform}/{protocol}/candidate_specific/alignment_{scatteritem}.bam",
        alignment_index="resources/{platform}/{protocol}/candidate_specific/alignment_{scatteritem}.bam.bai",
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
    shell:
        """
        cd {input.varlo}
        if [[ "{wildcards.platform}" == "Illumina_pe" || "{wildcards.platform}" == "Illumina_se" ]]; then
            PLATFORM="Illumina"
        else
            PLATFORM="{wildcards.platform}"
        fi
        echo $PLATFORM
        cargo run --release -- preprocess variants --omit-mapq-adjustment {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --read-type $PLATFORM --max-depth 5000 > {output}
        """
        # cargo run --release -- preprocess variants {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --read-type $PLATFORM --max-depth 5000 > {output}


rule call_methylation:
    input:
        varlo=directory("resources/tools/varlociraptor"),
        preprocess_obs="results/{platform}/{protocol}/normal_{scatteritem}.bcf",
        scenario="resources/scenario.yaml",
    output:
        "results/{platform}/{protocol}/calls_{scatteritem}.bcf",
    log:
        "logs/call_methylation_{platform}_{protocol}_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        """ 
        cd {input.varlo}
        cargo run --release -- call variants generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output}
        """


# TODO: Reactivate, right now it deletes too much data
rule filter_calls:
    input:
        varlo=directory("resources/tools/varlociraptor"),
        bcf="results/{platform}/{protocol}/calls_{scatteritem}.bcf",
    output:
        "results/{platform}/{protocol}/calls_{scatteritem}.filtered.bcf",
    log:
        "logs/filter_calls_{platform}_{protocol}_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        event="PRESENT",
    shell:
        """
        cd {input.varlo}
        cargo run --release -- filter-calls control-fdr --mode local-smart {input.bcf} --events {params.event} --fdr 0.005 > {output}
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
    shell:
        "cat {input} > {output}"


rule rename_varlo_output:
    input:
        "results/{platform}/{protocol}/varlo.vcf",
    output:
        "results/{platform}/{protocol}/result_files/varlo.bed",
    shell:
        "mv {input} {output}"
