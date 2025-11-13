rule compute_pandas_df:
    input:
        tool="results/{call_type}/{seq_platform}/called/{sample}/result_files/{method}.bed",
    output:
        "results/{call_type}/{seq_platform}/called/{sample}/result_files/{method}.parquet",
    conda:
        "../envs/plot.yaml"
    wildcard_constraints:
        method="(?!varlo|sample_df).*",
    log:
        "logs/plot_results/compute_pandas_df/{call_type}_{seq_platform}_{sample}_{method}.log",
    resources:
        mem_mb=64000,
    script:
        "../scripts/pandas_df_from_meth_output.py"


rule compute_varlo_df:
    input:
        tool="results/{call_type}/{seq_platform}/called/{sample}/result_files/varlo.bed",
    output:
        "results/{call_type}/{seq_platform}/{fdr}/{sample}/result_files/varlo.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/compute_varlo_df/{call_type}_{seq_platform}_{fdr}_{sample}.log",
    params:
        alpha=lambda wildcards: wildcards.fdr,
    resources:
        mem_mb=64000,
    script:
        "../scripts/pandas_df_from_meth_output.py"


# Computes one common df out of all single methylation caller dfs
rule common_tool_df:
    input:
        tools=lambda wildcards: expand(
            "results/{{call_type}}/{{seq_platform}}/called/{{sample}}/result_files/{method}.parquet",
            method=config["ref_tools"].get(wildcards.seq_platform, []),
        ),
        varlo="results/{call_type}/{seq_platform}/{fdr}/{sample}/result_files/varlo.parquet",
    output:
        sample_df="results/{call_type}/{seq_platform}/{fdr}/{sample}/result_files/sample_df.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/common_tool_df/{call_type}_{seq_platform}_{fdr}_{sample}.log",
    params:
        plot_type=config["plot_type"],
    wildcard_constraints:
        call_type="(?!multi_sample).*",
    resources:
        mem_mb=64000,
    script:
        "../scripts/common_tool_df.py"


# Plot overview tables over all samples of a sequencing platform with correlation information
# We correlate
# - all tools against each other on the same replicate
# - all replicates against each other for the same tool
rule compute_correlation:
    input:
        samples=lambda wildcards: expand(
            "results/{call_type}/{seq_platform}/{fdr}/{sample}/result_files/sample_df.parquet",
            call_type=wildcards.call_type,
            seq_platform=wildcards.seq_platform,
            fdr=wildcards.fdr,
            sample=config["data"].get(wildcards.seq_platform, wildcards.seq_platform),
        ),
    output:
        table="results/{call_type}/{seq_platform}/{fdr}/plots/replicates.hd5",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/compute_correlation/{call_type}_{seq_platform}_{fdr}.log",
    wildcard_constraints:
        call_type="(?!multi_sample).*",
    resources:
        mem_mb=64000,
    params:
        # plot_type=config["plot_type"],
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + ["varlo"],
        correlation_methods=config["correlation_methods"],
    script:
        "../scripts/compute_correlation.py"


rule replicates_heatmap:
    input:
        "results/{call_type}/{seq_platform}/{fdr}/plots/replicates.hd5",
    output:
        bias=report(
            "results/{call_type}/{seq_platform}/{fdr}/plots/{sample}_bias.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: f"{wildcards.seq_platform}",
            labels={
                "file": "bias",
                "sample": "{sample}",
                "fdr": "{fdr}",
            },
        ),
        heatmap=report(
            "results/{call_type}/{seq_platform}/{fdr}/plots/{sample}_heatmap.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: f"{wildcards.seq_platform}",
            labels={
                "file": "heatmap",
                "sample": "{sample}",
                "fdr": "{fdr}",
            },
        ),
    conda:
        "../envs/plot.yaml"
    resources:
        mem_mb=64000,
    log:
        "logs/plot_results/replicates_heatmap/{call_type}_{seq_platform}_{fdr}_{sample}_{plot_type}.log",
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + ["varlo"],
        sample=lambda wildcards: wildcards.sample,
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        fdr=lambda wildcards: wildcards.fdr,
        paper_plots=True,
    script:
        "../scripts/plot_heatmap_comparison.py"


# Compute common heatmap over all Illumina samples
rule heatmap_illumina_samples:
    input:
        "results/single_sample/Illumina_pe/{fdr}/plots/replicates.hd5",
    output:
        heatmap=report(
            "results/single_sample/Illumina_pe/{fdr}/plots/heatmap_all_samples.{plot_type}",
            category="single_sample",
            subcategory="Illumina_pe",
            labels={
                "file": "heatmap",
                "sample": "all samples",
                "fdr": "{fdr}",
            },
        ),
        bias=report(
            "results/single_sample/Illumina_pe/{fdr}/plots/bias_all_samples.{plot_type}",
            category="single_sample",
            subcategory="Illumina_pe",
            labels={
                "file": "bias",
                "sample": "all samples",
                "fdr": "{fdr}",
            },
        ),
        bar_plot_single_samples=report(
            "results/single_sample/Illumina_pe/{fdr}/plots/bar_plot_single_samples.{plot_type}",
            category="single_sample",
            subcategory="Illumina_pe",
            labels={
                "file": "bar_plot_single_samples",
                "sample": "all samples",
                "fdr": "{fdr}",
            },
        ),
    conda:
        "../envs/plot.yaml"
    resources:
        mem_mb=64000,
    log:
        "logs/plot_results/heatmap_illumina_samples/single_sample_Illumina_pe_{fdr}_{plot_type}.log",
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get("Illumina_pe", [])
        + ["varlo"],
        sample=lambda wildcards: config["samples"].get("Illumina_pe", []),
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        fdr=lambda wildcards: wildcards.fdr,
        paper_plots=True,
    script:
        "../scripts/plot_heatmap_comparison.py"


rule plot_runtime_comparison:
    input:
        benchmarks="benchmarks",
    output:
        tools="results/runtime_comparison_tools.{plot_type}",
        varlo="results/runtime_comparison_varlo.{plot_type}",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/plot_runtime_comparison/{plot_type}.log",
    script:
        "../scripts/plot_runtime_comparison.py"
