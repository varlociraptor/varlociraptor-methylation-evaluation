rule compute_pandas_df:
    input:
        tool="results/{call_type}/{seq_platform}/{protocol}/result_files/{method}.bed",
        # true_meth=lambda wildcards: expand(
        #     "resources/bed_avg_{chrom}.bedGraph",
        #     chrom=chromosome_by_seq_platform[wildcards.seq_platform],
        # ),
        # coverage="resources/{call_type}/{seq_platform}/{protocol}/cov.regions.bed",
    output:
        "results/{call_type}/{seq_platform}/{protocol}/result_files/{method}.parquet",
    conda:
        "../envs/plot.yaml"
    wildcard_constraints:
        protocol="(?!simulated_data).*",
        method="(?!varlo|protocol_df).*",
    log:
        "logs/plots/{call_type}/{seq_platform}/{protocol}/compute_pandas_df_{method}.log",
    params:
        # cov_bin_size=lambda wildcards: config["cov_bin_size"][wildcards.seq_platform],
        # cov_bins=lambda wildcards: config["cov_bins"][wildcards.seq_platform],
        # meth_type=config["meth_type"],
        # simulated=False,
        # prob_pres_threshhold=config["prob_pres_threshhold"],
        # prob_absent_threshhold=config["prob_absent_threshhold"],
    script:
        "../scripts/df_from_calls.py"


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
        # cov_bin_size=lambda wildcards: config["cov_bin_size"][wildcards.seq_platform],
        # cov_bins=lambda wildcards: config["cov_bins"][wildcards.seq_platform],
        # meth_type=config["meth_type"],
        # simulated=False,
        alpha=lambda wildcards: wildcards.fdr,
        # prob_pres_threshhold=config["prob_pres_threshhold"],
        # prob_absent_threshhold=config["prob_absent_threshhold"],
    script:
        "../scripts/df_from_calls.py"


# rule compute_pandas_df_simulated:
#     input:
#         tool="results/Illumina_pe/simulated_data/result_files/{method}.bed",
#         true_meth=expand(
#             "resources/Illumina_pe/simulated_data/chromosome_{chrom}_truth.bed",
#             chrom=config["seq_platforms"]["Illumina_pe"],
#         ),
#     output:
#         "results/Illumina_pe/simulated_data/result_files/{method}.parquet",
#     conda:
#         "../envs/plot.yaml"
#     params:
#         cov_bin_size=lambda wildcards: config["cov_bin_size"]["Illumina_pe"],
#         cov_bins=lambda wildcards: config["cov_bins"]["Illumina_pe"],
#         meth_type=config["meth_type"],
#         simulated=True,
#     log:
#         "logs/plots/Illumina_pe/simulated_data/compute_pandas_df_simulated_{method}.log",
#     script:
#         "../scripts/df_from_calls.py"


# rule plot_results_cov_specific:
#     input:
#         tool="results/{call_type}/{seq_platform}/{protocol}/result_files/{method}.parquet",
#     output:
#         plot=report(
#             "results/{call_type}/{seq_platform}/{protocol}/plots/{method}/plots_{cov_bin, [0-9]+}.{plot_type}",
#             category=lambda wildcards: f"{wildcards.seq_platform} - {wildcards.protocol}",
#             subcategory=lambda wildcards: f"{wildcards.method} plots",
#             labels=lambda wildcards: {
#                 "coverage": f"{int(wildcards.cov_bin)* int(config['cov_bin_size'][wildcards.seq_platform])} - {int(wildcards.cov_bin)* int(config['cov_bin_size'][wildcards.seq_platform])+ int(config['cov_bin_size'][wildcards.seq_platform])}",
#             },
#         ),
#     conda:
#         "../envs/plot.yaml"
#     log:
#         "logs/plots/{call_type}/{seq_platform}/{protocol}/results_cov_specific_{method}_{cov_bin}_{plot_type}.log",
#     params:
#         plot_type=config["plot_type"],
#         cov_bin=lambda wildcards: wildcards.cov_bin,
#     script:
#         "../scripts/scatter_plot.py"


# rule plot_results_all_cov:
#     input:
#         tool="results/{call_type}/{seq_platform}/{protocol}/result_files/{method}.parquet",
#     output:
#         plot=report(
#             "results/{call_type}/{seq_platform}/{protocol}/plots/{method}/plots_all.{plot_type}",
#             category=lambda wildcards: f"{wildcards.seq_platform} - {wildcards.protocol}",
#             subcategory=lambda wildcards: f"{wildcards.method} plots",
#             labels=lambda wildcards: {
#                 "coverage": "all",
#             },
#         ),
#     conda:
#         "../envs/plot.yaml"
#     log:
#         "logs/plots/{call_type}/{seq_platform}/{protocol}/all_cov_{method}_{plot_type}.log",
#     params:
#         plot_type=config["plot_type"],
#         cov_bin=-1,
#     script:
#         "../scripts/scatter_plot.py"


# # TODO: Catch JavaScript out of Memory exception and give advise to use html for bigger datasets
# rule scatter_ref_method:
#     input:
#         varlo="results/{call_type}/{seq_platform}/{protocol}/result_files/varlo.parquet",
#         ref_tool="results/{call_type}/{seq_platform}/{protocol}/result_files/{method}.parquet",
#     output:
#         scatter_plot=report(
#             expand(
#                 "results/{{seq_platform}}/{{protocol}}/plots/scatter_comp_{{method}}.{plot_type}",
#                 plot_type=config["plot_type"],
#             ),
#             category=lambda wildcards: f"{wildcards.seq_platform} - {wildcards.protocol}",
#             subcategory="All Comparisions",
#             labels=lambda wildcards: {
#                 "type": "scatter comparision",
#                 "method": wildcards.method,
#             },
#         ),
#     conda:
#         "../envs/plot.yaml"
#     log:
#         "logs/plots/{call_type}/{seq_platform}/{protocol}/tool_comparision_{method}.log",
#     params:
#         plot_type=config["plot_type"],
#         cov_bins=config["cov_bins"],
#     resources:
#         mem_mb=8000
#     script:
#         "../scripts/plot_scatter_compare.py"


# rule compute_precision_recall:
#     input:
#         tool="results/{call_type}/{seq_platform}/{protocol}/result_files/{method}.parquet",
#     output:
#         precall="results/{call_type}/{seq_platform}/{protocol}/plots/{method}/precall_{cov_bin}.csv",
#     conda:
#         "../envs/plot.yaml"
#     params:
#         plot_type=config["plot_type"],
#         cov_bin=lambda wildcards: wildcards.cov_bin,
#         cov_bin_size=lambda wildcards: config["cov_bin_size"][wildcards.seq_platform],
#         meth_type=config["meth_type"],
#         meth_threshold=config["meth_threshold"],
#     log:
#         "logs/plots/{call_type}/{seq_platform}/{protocol}/compute_precision_recall_{method}_{cov_bin}.log",
#     script:
#         "../scripts/compute_precision_recall.py"


# rule plot_precision_recall:
#     input:
#         get_precision_recall_csvs,
#     output:
#         report(
#             expand(
#                 "results/{{seq_platform}}/{{protocol}}/plots/precall.{{plot_type}}",
#             ),
#             category=lambda wildcards: f"{wildcards.seq_platform} - {wildcards.protocol}",
#             subcategory="All Comparisions",
#             labels=lambda wildcards: {"type": "precision recall", "method": "all"},
#         ),
#     conda:
#         "../envs/plot.yaml"
#     log:
#         "logs/plots/{call_type}/{seq_platform}/{protocol}/plot_precision_recall_{plot_type}.log",
#     params:
#         cov_bin=lambda wildcards: config["cov_bins"][wildcards.seq_platform],
#         plot_type=config["plot_type"],
#     script:
#         "../scripts/plot_precision_recall.py"


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
rule compute_correlation_tables:
    input:
        samples=lambda wildcards: expand(
            "results/{call_type}/{seq_platform}/{fdr}/{protocol}/result_files/protocol_df.parquet",
            call_type=wildcards.call_type,
            seq_platform=wildcards.seq_platform,
            fdr=wildcards.fdr,
            protocol=config["data"][wildcards.seq_platform],
        ),
    output:
        # plot=report(
        #     "results/{call_type}/{seq_platform}/{fdr}/plots/correlation_table.{plot_type}",
        #     category="{seq_platform}",
        #     subcategory=lambda wildcards: f"{wildcards.fdr}",
        #     labels={"file": "correlation comparision"},
        # ),
        table="results/{call_type}/{seq_platform}/{fdr}/plots/replicates.hd5",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/{call_type}/{seq_platform}/{fdr}/correlation_tables.log",
    wildcard_constraints:
        call_type="(?!multi_sample).*",
    params:
        # plot_type=config["plot_type"],
        meth_callers=lambda wildcards: config["ref_tools"].get(wildcards.seq_platform, []) + ["varlo"],
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
        benchmarks="benchmarks"
    output:
        tools="results/runtime_comparison_tools.{plot_type}",
        varlo="results/runtime_comparison_varlo.{plot_type}",
        # memory="results/memory_comparison.html"
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/plot_runtime_comparisons.py"
