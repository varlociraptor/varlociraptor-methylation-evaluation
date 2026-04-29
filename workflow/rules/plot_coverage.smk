# rule filter_bam_mapq:
#     input:
#         bam="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
#         bai="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
#     output:
#         bam="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.mapq_{mapq}.bam",
#         bai="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.mapq_{mapq}.bam.bai",
#     threads: 4
#     log:
#         "logs/mapq60/{seq_platform}_{sample}_{mapq}.log"
#     params:
#         mapq=lambda wildcards: f" -q {wildcards.mapq}" if wildcards.mapq != "all" else ""
#     shell:
#         """
#         samtools view -@ {threads} -b {params.mapq} {input.bam} > {output.bam} 2> {log}
#         samtools index -@ {threads} {output.bam}
#         """

rule compute_coverage:
    input:
        bam="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        bai="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        bed=lambda wildcards: expand(
            "resources/{chrom}/candidates.bed",
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
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v5.5.2/bio/mosdepth"

# rule unzip_coverage:
#     input:
#         "results/{call_type}/{seq_platform}/coverages/{sample}.regions.bed{mapq}.gz",
#     output:
#         "results/{call_type}/{seq_platform}/coverages/{sample}.regions{mapq}.bed",
#     log:
#         "logs/mason/mason_unzip_coverage/{call_type}_{seq_platform}_{sample}_{mapq}.log",
#     conda:
#         "../envs/general.yaml"
#     shell:
#         "gunzip -c {input} > {output} 2> {log}"



rule coverage_plots:
    input:
        coverage="results/{call_type}/{seq_platform}/coverages/{sample}.regions.bed.gz",
        meth_data="results/{call_type}/{seq_platform}/result_files/sample_df_{sample}.parquet",
    output:
        meth_level_to_cov="results/{call_type}/{seq_platform}/plots/{sample}_meth_level_to_cov.{plot_type}",
        coverage_retained=report(
            "results/{call_type}/{seq_platform}/plots/{sample}_coverage_retained.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: f"{wildcards.seq_platform}",
            labels={
                "file": "coverage_retained",
                "sample": "{sample}",
            },
            caption="../report/coverage_retained.rst",
        ),
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=4000,
    log:
        "logs/plot_results/coverage_plots/{call_type}_{seq_platform}_{sample}_{plot_type}.log",
    params:
        sample=lambda wildcards: config["samples"].get(wildcards.seq_platform, []),
        plot_type=lambda wildcards: wildcards.plot_type,
    script:
        "../scripts/plot_coverage_retained.py"

rule stratify_mae:
    input:
        coverage_all_01="results/{call_type}/{seq_platform}/coverages/{sample}_REP01_all.regions.bed.gz",
        coverage_all_02="results/{call_type}/{seq_platform}/coverages/{sample}_REP02_all.regions.bed.gz",
        coverage_60_01="results/{call_type}/{seq_platform}/coverages/{sample}_REP01_60.regions.bed.gz",
        coverage_60_02="results/{call_type}/{seq_platform}/coverages/{sample}_REP02_60.regions.bed.gz",
        # coverage_60="results/{call_type}/{seq_platform}/coverages/{sample}.regions.bed.mapq60.gz",
        meth_data="results/{call_type}/{seq_platform}/result_files/replicates.parquet",
    output:
        mae="results/{call_type}/{seq_platform}/plots/{sample}_mae.{plot_type}",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=4000,
    log:
        "logs/plot_results/stratify_mae/{call_type}_{seq_platform}_{sample}_{plot_type}.log",
    params:
        sample=lambda wildcards: wildcards.sample,
        plot_type=lambda wildcards: wildcards.plot_type,
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
    script:
        "../scripts/plot_stratify_mae.py"
