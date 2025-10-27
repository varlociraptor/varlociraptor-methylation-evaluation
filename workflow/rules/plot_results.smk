rule compute_pandas_df:
    input:
        tool="results/{call_type}/{seq_platform}/called/{sample}/result_files/{method}.bed",
    output:
        "results/{call_type}/{seq_platform}/called/{sample}/result_files/{method}.parquet",
    conda:
        "../envs/plot.yaml"
    wildcard_constraints:
        sample="(?!simulated_data).*",
        method="(?!varlo|sample_df).*",
    log:
        "logs/plots/{call_type}/{seq_platform}/{sample}/compute_pandas_df_{method}.log",
    script:
        "../scripts/pandas_df_from_meth_output.py"


rule compute_varlo_df:
    input:
        tool="results/{call_type}/{seq_platform}/called/{sample}/result_files/varlo.bed",
    output:
        "results/{call_type}/{seq_platform}/{fdr}/{sample}/result_files/varlo.parquet",
    conda:
        "../envs/plot.yaml"
    wildcard_constraints:
        sample="(?!simulated_data).*",
    log:
        "logs/plots/{call_type}/{seq_platform}/{fdr}/{sample}/compute_pvarlo_df.log",
    params:
        alpha=lambda wildcards: wildcards.fdr,
    script:
        "../scripts/pandas_df_from_meth_output.py"


# Computes one common df out of all single method dfs
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
        "logs/plots/{call_type}/{seq_platform}/{fdr}/{sample}/common_tool_df.log",
    params:
        plot_type=config["plot_type"],
    wildcard_constraints:
        call_type="(?!multi_sample).*",
    resources:
        mem_mb=16000,
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
        "logs/plots/{call_type}/{seq_platform}/{fdr}/correlation_tables.log",
    wildcard_constraints:
        call_type="(?!multi_sample).*",
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
        report(
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
        mem_mb=32000,
    log:
        "logs/plots/{call_type}/{seq_platform}/{fdr}/plot_heatmap_comparision_{sample}_{plot_type}.log",
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + ["varlo"],
        # meth_callers=["varlo"],
        sample=lambda wildcards: wildcards.sample,
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        # correlation_method=config["correlation_method"],
    script:
        "../scripts/plot_heatmap_comparision.py"


# Compute common heatmap over all Illumina samples
rule heatmap_illumina_samples:
    input:
        "results/{call_type}/{seq_platform}/{fdr}/plots/replicates.hd5",
    output:
        report(
            "results/{call_type}/{seq_platform}/{fdr}/plots/heatmap_all_samples.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: f"{wildcards.seq_platform}",
            labels={
                "file": "heatmap",
                "sample": "all samples",
                "fdr": "{fdr}",
            },
        ),
    conda:
        "../envs/plot.yaml"
    resources:
        mem_mb=32000,
    log:
        "logs/plots/{call_type}/{seq_platform}/{fdr}/plot_heatmap_comparision_{plot_type}.log",
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + ["varlo"],
        # sample="all",
        sample=lambda wildcards: config["samples"].get(wildcards.seq_platform, []),
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        # correlation_method=config["correlation_method"],
    script:
        "../scripts/plot_heatmap_comparision.py"


rule plot_runtime_comparison:
    input:
        benchmarks="benchmarks",
    output:
        tools="results/runtime_comparison_tools.{plot_type}",
        varlo="results/runtime_comparison_varlo.{plot_type}",
        # memory="results/memory_comparison.html"
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/plot_runtime_comparision_{plot_type}.log",
    script:
        "../scripts/plot_runtime_comparision.py"
