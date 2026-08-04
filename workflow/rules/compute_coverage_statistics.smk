# Our candidate file spans CG positions, we are only interested on the coverage of C positions
rule focus_candidates_on_c:
    input:
        "resources/{chrom}/candidates.bed",
    output:
        "resources/{chrom}/candidates_focus_c.bed",
    log:
        "logs/mason/focus_candidates_on_c/{chrom}.log",
    conda:
        "../envs/general.yaml"
    shell:
        r"""
        awk 'BEGIN{{OFS="\t"}}{{
            center=$2+1
            print $1, center, center+1
        }}' {input} > {output}
        """

rule compute_coverage:
    input:
        bam="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        bai="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        bed=lambda wildcards: expand(
            "resources/{chrom}/candidates_focus_c.bed",
            chrom=config["seq_platforms"].get(wildcards.seq_platform, []),
        ),
    output:
        "results/{call_type}/{seq_platform}/coverages/{sample}_{mapq}.mosdepth.global.dist.txt",
        "results/{call_type}/{seq_platform}/coverages/{sample}_{mapq}.mosdepth.region.dist.txt",
        "results/{call_type}/{seq_platform}/coverages/{sample}_{mapq}.regions.bed.gz",
        summary="results/{call_type}/{seq_platform}/coverages/{sample}_{mapq}.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mason/mason_coverage/{call_type}_{seq_platform}_{sample}_{mapq}.log",
    params:
        extra=lambda wildcards: "--no-per-base --use-median" + f" --mapq {wildcards.mapq}" if wildcards.mapq != "all" else "",  # optional
    wildcard_constraints:
        sample="(?!all_samples).*",
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v5.5.2/bio/mosdepth"

rule common_coverage_Illumina:
    input:
        lambda wildcards: expand("results/single_sample/Illumina_pe/coverages/{sample}_{{REP}}_{{mapq}}.regions.bed.gz", sample=config["samples"]["Illumina_pe"])
    output:
        "results/single_sample/Illumina_pe/coverages/all_samples_{REP}_{mapq}.regions.bed.gz",
    threads: 1
    conda:
        "../envs/python.yaml"
    log:
        "logs/common_coverage_Illumina/all_samples_{REP}_{mapq}.log",
    script:
        "../scripts/merge_coverage.py"


rule compute_coverage_retention:
    input:
        coverage_01="results/{call_type}/{seq_platform}/coverages/{sample}_REP01_{mapq}.regions.bed.gz",
        coverage_02="results/{call_type}/{seq_platform}/coverages/{sample}_REP02_{mapq}.regions.bed.gz",
        meth_data="results/{call_type}/{seq_platform}/result_files/replicates.parquet",
    output:
        mae=report("results/{call_type}/{seq_platform}/plots/{sample}_dist_{mapq}.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: seq_platform_to_name[wildcards.seq_platform],
            labels= {
                "file_type": "coverage_retention",
                "sample": "{sample}",
            },
            caption="../report/coverage_retained.rst",
        ),
        parquet="results/{call_type}/{seq_platform}/coverages/{sample}_coverage_plot_{plot_type}_{mapq}.parquet",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=16000,
    log:
        "logs/plot_results/stratify_mae/{call_type}_{seq_platform}_{sample}_{mapq}_{plot_type}.log",
    params:
        sample=lambda wildcards: wildcards.sample,
        plot_type=lambda wildcards: wildcards.plot_type,
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + [f"varlo_{fdr}" for fdr in config["fdr_alpha"]],
        # How many datapoints to include in the plot (quantile). If 100 we have really high coverages and can't see shit
        quantile=lambda wildcards: 0.995 if wildcards.seq_platform == "Illumina_pe" else 0.95,
        bin_size=5,
    script:
        "../scripts/compute_coverage_retention.py"
