# TODO: Missingoutputexception
rule download_varlociraptor:
    output:
        "resources/tools/varlociraptor/Cargo.toml",
    log:
        "logs/varlociraptor/download.log",
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
        "resources/tools/varlociraptor",
    output:
        "resources/tools/varlociraptor/target/debug/varlociraptor",
    log:
        "logs/varlociraptor/build.log",
    conda:
        "../envs/varlociraptor.yaml"
    resources:
        mem_mb=8000,
    shell:
        """
        echo {output} >> {log}
        cd {input}
        CARGO_TARGET_DIR={resources.tmpdir}/cargo_target \
            cargo build >> {log}
        mkdir -p $(dirname {output})
        cp {resources.tmpdir}/cargo_target/debug/varlociraptor {output}
        echo "Show output file:" >> {log}
        ls -ahl $(dirname {output}) >> {log}
        """
        # mkdir -p $(dirname {output})
        # touch {output}
        # ls $(dirname {output}) >> {log}


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
            "resources/{chrom}/candidates_{{scatteritem}}.bcf",
            chrom=chromosome_by_seq_platform[wildcards.seq_platform],
        ),
    output:
        "results/{seq_platform}/{protocol}/normal_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/{seq_platform}/{protocol}/compute_meth_observations_{scatteritem}.log",
    benchmark:
        "benchmarks/{seq_platform}/preprocessing/{protocol}_{scatteritem}.bwa.benchmark.txt"
    conda:
        "../envs/varlociraptor.yaml"
    resources:
        mem_mb=16000,
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
        "results/single_sample/{seq_platform}/{protocol}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/{seq_platform}/{protocol}/call_methylation_{scatteritem}.log",
    benchmark:
        "benchmarks/{seq_platform}/calling/single_sample/{protocol}_{scatteritem}.bwa.benchmark.txt"
    conda:
        "../envs/varlociraptor.yaml"
    wildcard_constraints:
        seq_platform="(?!common_calls).*",
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output} 2> {log}"


# TODO: Reactivate, right now it deletes too much data
# rule filter_calls:
#     input:
#         varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
#         bcf="results/{seq_platform}/{protocol}/calls_{scatteritem}.bcf",
#     output:
#         "results/{seq_platform}/{fdr}/{protocol}/calls_{scatteritem}.bcf",
#     log:
#         "logs/varlociraptor/{seq_platform}/{fdr}/{protocol}/filter_calls_{scatteritem}.log",
#     conda:
#         "../envs/varlociraptor.yaml"
#     params:
#         event="LOW",
#         fdr_alpha=config["fdr_alpha"]
#     shell:
#         "{input.varlo} filter-calls control-fdr --mode local-smart {input.bcf} --events LOW HIGH --fdr {params.fdr_alpha} > {output} 2> {log}"


# TODO: Skip this step, right now it would be useless since I debug so much
rule calls_to_vcf:
    input:
        "results/{call_type}/{seq_platform}/{protocol}/calls_{scatteritem}.bcf",
    output:
        "results/{call_type}/{seq_platform}/{protocol}/calls_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/varlociraptor/{call_type}/{seq_platform}/{protocol}/calls_to_vcf_{scatteritem}.log",
    threads: 10
    shell:
        # "touch {output} 2> {log}"
        "bcftools view --threads {threads} {input} -o {output} 2> {log}"


rule gather_calls:
    input:
        gather.split_candidates(
            "results/{{call_type}}/{{seq_platform}}/{{protocol}}/calls_{scatteritem}.vcf"
        ),
    output:
        "results/{call_type}/{seq_platform}/{protocol}/result_files/varlo.bed",
    log:
        "logs/varlociraptor/{call_type}/{seq_platform}/{protocol}/gather_calls.log",
    shell:
        "cat {input} > {output} 2> {log}"


# rule rename_varlo_output:
#     input:
#         "results/{seq_platform}/{protocol}/varlo.vcf",
#     output:
#         "results/{seq_platform}/{protocol}/result_files/varlo.bed",
#     log:
#         "logs/varlociraptor/{seq_platform}/{protocol}/rename_varlo_output.log",
#     shell:
#         "mv {input} {output}"
