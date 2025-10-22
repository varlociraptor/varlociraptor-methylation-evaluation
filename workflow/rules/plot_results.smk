rule compute_pandas_df:
    input:
        tool="results/{call_type}/{seq_platform}/{protocol}/result_files/{method}.bed",
    output:
        "results/{call_type}/{seq_platform}/{protocol}/result_files/{method}.parquet",
    conda:
        "../envs/plot.yaml"
    wildcard_constraints:
        protocol="(?!simulated_data).*",
        method="(?!varlo|protocol_df).*",
    log:
        "logs/plots/{call_type}/{seq_platform}/{protocol}/compute_pandas_df_{method}.log",
    script:
        "../scripts/pandas_df_from_meth_output.py"


rule compute_varlo_df:
    input:
        tool="results/{call_type}/{seq_platform}/{protocol}/result_files/varlo.bed",
    output:
        "results/{call_type}/{seq_platform}/{fdr}/{protocol}/result_files/varlo.parquet",
    conda:
        "../envs/plot.yaml"
    wildcard_constraints:
        protocol="(?!simulated_data).*",
    log:
        "logs/plots/{call_type}/{seq_platform}/{fdr}/{protocol}/compute_pvarlo_df.log",
    params:
        alpha=lambda wildcards: wildcards.fdr,
    script:
        "../scripts/pandas_df_from_meth_output.py"


# Computes one common df out of all single method dfs
rule common_tool_df:
    input:
        tools=lambda wildcards: expand(
            "results/{{call_type}}/{{seq_platform}}/{{protocol}}/result_files/{method}.parquet",
            method=config["ref_tools"][wildcards.seq_platform],
        ),
        varlo="results/{call_type}/{seq_platform}/{fdr}/{protocol}/result_files/varlo.parquet",
    output:
        protocol_df="results/{call_type}/{seq_platform}/{fdr}/{protocol}/result_files/protocol_df.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/{call_type}/{seq_platform}/{fdr}/{protocol}/common_tool_df.log",
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
            "results/{call_type}/{seq_platform}/{fdr}/{protocol}/result_files/protocol_df.parquet",
            call_type=wildcards.call_type,
            seq_platform=wildcards.seq_platform,
            fdr=wildcards.fdr,
            protocol=config["data"][wildcards.seq_platform],
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
        "../scripts/compute_correlation_tables.py"


rule replicates_heatmap:
    input:
        "results/{call_type}/{seq_platform}/{fdr}/plots/replicates.hd5",
    output:
        report(
            "results/{call_type}/{seq_platform}/{fdr}/plots/{protocol}_heatmap.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: f"{wildcards.seq_platform}",
            labels={
                "file": "heatmap",
                "protocol": "{protocol}",
                "fdr": "{fdr}",
            },
        ),
    conda:
        "../envs/plot.yaml"
    resources:
        mem_mb=32000,
    log:
        "logs/plots/{call_type}/{seq_platform}/{fdr}/heatmap_replicates_{protocol}_{plot_type}.log",
    params:
        # meth_callers = lambda wildcards: config["ref_tools"].get(wildcards.seq_platform, []) + ["varlo"],
        meth_callers=["varlo"],
        protocol=lambda wildcards: wildcards.protocol,
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        # correlation_method=config["correlation_method"],
    script:
        "../scripts/heatmap_replicates.py"


# Compute common heatmap over all Illumina protocols
rule heatmap_illumina_protocols:
    input:
        "results/{call_type}/{seq_platform}/{fdr}/plots/replicates.hd5",
    output:
        report(
            "results/{call_type}/{seq_platform}/{fdr}/plots/heatmap_all_protocols.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: f"{wildcards.seq_platform}",
            labels={
                "file": "heatmap",
                "protocol": "all protocols",
                "fdr": "{fdr}",
            },
        ),
    conda:
        "../envs/plot.yaml"
    resources:
        mem_mb=32000,
    log:
        "logs/plots/{call_type}/{seq_platform}/{fdr}/heatmap_replicates_{plot_type}.log",
    params:
        meth_callers=lambda wildcards: config["ref_tools"][wildcards.seq_platform]
        + ["varlo"],
        # protocol="all",
        protocol=lambda wildcards: config["protocols"][wildcards.seq_platform],
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        # correlation_method=config["correlation_method"],
    script:
        "../scripts/heatmap_replicates.py"


rule plot_runtime_comparison:
    input:
        benchmarks="benchmarks",
    output:
        tools="results/runtime_comparison_tools.{plot_type}",
        varlo="results/runtime_comparison_varlo.{plot_type}",
        # memory="results/memory_comparison.html"
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/plot_runtime_comparisons.py"
