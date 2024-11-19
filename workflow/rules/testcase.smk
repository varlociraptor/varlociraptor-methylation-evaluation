rule preprocess_testcase:
    input:
        chromosome=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
        genome_index=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
        alignments="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        candidates=lambda wildcards: expand(
            "resources/{chrom}/candidates.bcf",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
    output:
        # True output is in varlociraptor 
        temp("results/{platform}/{protocol}/preprocess_test_cration_finished.txt"),
    log:
        "logs/preprocess_testcase_{platform}_{protocol}.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        varlo_path=config["varlo_path"],
        pipeline_path=config["pipeline_path"],
        locus=config["testcase"]["locus"],
        prefix=config["testcase"]["prefix"],
    shell:
        """ 
        cd {params.varlo_path}
        cargo run --release -- preprocess variants --testcase-locus {params.locus} --testcase-prefix {params.prefix} {params.pipeline_path}{input.chromosome} --candidates {params.pipeline_path}{input.candidates} --bam {params.pipeline_path}{input.alignments} --read-type Illumina --max-depth 5000 > {params.pipeline_path}{output}
        touch {params.pipeline_path}{output}
        """


rule call_testcase:
    input:
        preprocess_obs="results/{platform}/{protocol}/normal_6-of-20.bcf",
        scenario="resources/scenario.yaml",
    output:
        # True output is in varlociraptor 
        temp("results/{platform}/{protocol}/call_test_creation_finished.txt"),
    log:
        "logs/call_testcase_{platform}_{protocol}.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        varlo_path=config["varlo_path"],
        pipeline_path=config["pipeline_path"],
        locus=config["testcase"]["locus"],
        prefix=config["testcase"]["prefix"],
    shell:
        """ 
        cd {params.varlo_path}
        cargo run --release -- call variants --testcase-locus {params.locus} --testcase-prefix {params.prefix} --omit-strand-bias generic --scenario {params.pipeline_path}{input.scenario} --obs normal={params.pipeline_path}{input.preprocess_obs}
        touch {params.pipeline_path}{output}
        """
