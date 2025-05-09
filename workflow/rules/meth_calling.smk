# TODO: Missingoutputexception
rule download_varlociraptor:
    output:
        "resources/tools/varlociraptor/Cargo.toml",
    log:
        "../logs/varlociraptor/download.log",
    conda:
        "../envs/shell_cmds.yaml"
    shell:
        """
        mkdir -p resources/tools 2> {log}
        cd resources/tools 2> {log}
        git clone https://github.com/varlociraptor/varlociraptor.git 2> {log}
        cd varlociraptor 2> {log}
        git checkout methylation-paired-end-master-new 2> {log}
        """


rule build_varlociraptor:
    input:
        "resources/tools/varlociraptor/Cargo.toml",
    output:
        "resources/tools/varlociraptor/target/debug/varlociraptor",
    log:
        "../logs/varlociraptor/build.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        """
        cd $(dirname {input})
        cargo build 2> {log}
        """


rule compute_meth_observations:
    input:
        varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
        chromosome=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=chromosome_by_seq_platform[wildcards.seq_platform],
        ),
        genome_index=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=chromosome_by_seq_platform[wildcards.seq_platform],
        ),
        alignments="resources/{seq_platform}/{protocol}/candidate_specific/alignment_{scatteritem}.bam",
        alignment_index="resources/{seq_platform}/{protocol}/candidate_specific/alignment_{scatteritem}.bam.bai",
        candidates=lambda wildcards: expand(
            "resources/{{seq_platform}}/{{protocol}}/{chrom}/candidates_{{scatteritem}}.bcf",
            chrom=chromosome_by_seq_platform[wildcards.seq_platform],
        ),
    output:
        "results/{seq_platform}/{protocol}/normal_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/{seq_platform}/{protocol}/compute_meth_observations_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        """
        if [[ "{wildcards.seq_platform}" == "Illumina_pe" || "{wildcards.seq_platform}" == "Illumina_se" ]]; then
            seq_platform="Illumina"
        else
            seq_platform="{wildcards.seq_platform}"
        fi
        {input.varlo} preprocess variants --omit-mapq-adjustment {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --read-type $seq_platform --max-depth 5000 > {output} 2> {log}
        """


rule call_methylation:
    input:
        varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
        preprocess_obs="results/{seq_platform}/{protocol}/normal_{scatteritem}.bcf",
        scenario="resources/scenario.yaml",
    output:
        "results/{seq_platform}/{protocol}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/{seq_platform}/{protocol}/call_methylation_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output} 2> {log}"


# TODO: Reactivate, right now it deletes too much data
rule filter_calls:
    input:
        varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
        bcf="results/{seq_platform}/{protocol}/calls_{scatteritem}.bcf",
    output:
        "results/{seq_platform}/{protocol}/calls_{scatteritem}.filtered.bcf",
    log:
        "logs/varlociraptor/{seq_platform}/{protocol}/filter_calls_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        event="PRESENT",
    shell:
        "{input.varlo} filter-calls control-fdr --mode local-smart {input.bcf} --events {params.event} --fdr 0.005 > {output} 2> {log}"


# TODO: Skip this step, right now it would be useless since I debug so much
rule calls_to_vcf:
    input:
        "results/{seq_platform}/{protocol}/calls_{scatteritem}.bcf",
    output:
        "results/{seq_platform}/{protocol}/calls_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/varlociraptor/{seq_platform}/{protocol}/calls_to_vcf_{scatteritem}.log",
    threads: 10
    shell:
        "bcftools view --threads {threads} {input} -o {output} 2> {log}"


rule gather_calls:
    input:
        gather.split_candidates(
            "results/{{seq_platform}}/{{protocol}}/calls_{scatteritem}.vcf"
        ),
    output:
        "results/{seq_platform}/{protocol}/varlo.vcf",
    log:
        "logs/varlociraptor/{seq_platform}/{protocol}/gather_calls.log",
    shell:
        "cat {input} > {output} 2> {log}"


rule rename_varlo_output:
    input:
        "results/{seq_platform}/{protocol}/varlo.vcf",
    output:
        "results/{seq_platform}/{protocol}/result_files/varlo.bed",
    log:
        "logs/varlociraptor/{seq_platform}/{protocol}/rename_varlo_output.log",
    shell:
        "mv {input} {output} 2> {log}"
