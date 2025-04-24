rule compute_pandas_df:
    input:
        tool="results/{platform}/{protocol}/result_files/{method}.bed",
        true_meth=lambda wildcards: expand(
            "resources/bed_avg_{chrom}.bedGraph",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
        coverage="resources/{platform}/{protocol}/cov.regions.bed",
    output:
        "results/{platform}/{protocol}/result_files/{method}.parquet",
    conda:
        "../envs/plot.yaml"
    wildcard_constraints:
        protocol="(?!simulated_data).*",
    log:
        "logs/plots/{platform}/{protocol}/compute_pandas_df_{method}.log",
    params:
        cov_bin_size=lambda wildcards: config["cov_bin_size"][wildcards.platform],
        cov_bins=lambda wildcards: config["cov_bins"][wildcards.platform],
        meth_type=config["meth_type"],
        simulated=False,
    script:
        "../scripts/df_from_calls.py"


rule compute_pandas_df_simulated:
    input:
        tool="results/Illumina_pe/simulated_data/result_files/{method}.bed",
        true_meth=expand(
            "resources/Illumina_pe/simulated_data/chromosome_{chrom}_truth.bed",
            chrom=config["platforms"]["Illumina_pe"],
        ),
    output:
        "results/Illumina_pe/simulated_data/result_files/{method}.parquet",
    conda:
        "../envs/plot.yaml"
    params:
        cov_bin_size=lambda wildcards: config["cov_bin_size"]["Illumina_pe"],
        cov_bins=lambda wildcards: config["cov_bins"]["Illumina_pe"],
        meth_type=config["meth_type"],
        simulated=True,
    log:
        "logs/plots/Illumina_pe/simulated_data/compute_pandas_df_simulated_{method}.log",
    script:
        "../scripts/df_from_calls.py"


rule plot_results_cov_specific:
    input:
        tool="results/{platform}/{protocol}/result_files/{method}.parquet",
    output:
        plot=report(
            "results/{platform}/{protocol}/plots/{method}/plots_{cov_bin, [0-9]+}.{plot_type}",
            category=lambda wildcards: f"{wildcards.platform} - {wildcards.protocol}",
            subcategory=lambda wildcards: f"{wildcards.method} plots",
            labels=lambda wildcards: {
                "coverage": f"{int(wildcards.cov_bin) * int(config['cov_bin_size'][wildcards.platform])} - {int(wildcards.cov_bin) * int(config['cov_bin_size'][wildcards.platform]) + int(config['cov_bin_size'][wildcards.platform])}",
            },
        ),
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/{platform}/{protocol}/results_cov_specific_{method}_{cov_bin}_{plot_type}.log",
    params:
        plot_type=config["plot_type"],
        prob_pres_threshhold=config["prob_pres_threshhold"],
        cov_bin=lambda wildcards: wildcards.cov_bin,
    script:
        "../scripts/scatter_plot.py"


rule plot_results_all_cov:
    input:
        tool="results/{platform}/{protocol}/result_files/{method}.parquet",
    output:
        plot=report(
            "results/{platform}/{protocol}/plots/{method}/plots_all.{plot_type}",
            category=lambda wildcards: f"{wildcards.platform} - {wildcards.protocol}",
            subcategory=lambda wildcards: f"{wildcards.method} plots",
            labels=lambda wildcards: {
                "coverage": "all",
            },
        ),
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/{platform}/{protocol}/all_cov_{method}_{plot_type}.log",
    params:
        plot_type=config["plot_type"],
        prob_pres_threshhold=config["prob_pres_threshhold"],
        cov_bin=-1,
    script:
        "../scripts/scatter_plot.py"


rule plot_scatter_comparision:
    input:
        varlo="results/{platform}/{protocol}/result_files/varlo.parquet",
        ref_tool="results/{platform}/{protocol}/result_files/{method}.parquet",
    output:
        scatter_plot=report(
            expand(
                "results/{{platform}}/{{protocol}}/plots/scatter_comp_{{method}}.{plot_type}",
                plot_type=config["plot_type"],
            ),
            category=lambda wildcards: f"{wildcards.platform} - {wildcards.protocol}",
            subcategory="All Comparisions",
            labels=lambda wildcards: {
                "type": "scatter comparision",
                "method": wildcards.method,
            },
        ),
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/{platform}/{protocol}/tool_comparision_{method}.log",
    params:
        plot_type=config["plot_type"],
        cov_bins=config["cov_bins"],
        prob_pres_threshhold=config["prob_pres_threshhold"],
    script:
        "../scripts/plot_scatter_compare.py"


rule compute_precision_recall:
    input:
        tool="results/{platform}/{protocol}/result_files/{method}.parquet",
    output:
        precall="results/{platform}/{protocol}/plots/{method}/precall_{cov_bin}.csv",
    conda:
        "../envs/plot.yaml"
    params:
        plot_type=config["plot_type"],
        cov_bin=lambda wildcards: wildcards.cov_bin,
        cov_bin_size=lambda wildcards: config["cov_bin_size"][wildcards.platform],
        prob_pres_threshhold=config["prob_pres_threshhold"],
        meth_type=config["meth_type"],
    log:
        "logs/plots/{platform}/{protocol}/compute_precision_recall_{method}_{cov_bin}.log",
    script:
        "../scripts/compute_precision_recall.py"


rule plot_precision_recall:
    input:
        get_precision_recall_csvs,
    output:
        report(
            expand(
                "results/{{platform}}/{{protocol}}/plots/precall.{{plot_type}}",
            ),
            category=lambda wildcards: f"{wildcards.platform} - {wildcards.protocol}",
            subcategory="All Comparisions",
            labels=lambda wildcards: {"type": "precision recall", "method": "all"},
        ),
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/{platform}/{protocol}/plot_precision_recall_{plot_type}.log",
    params:
        cov_bin=lambda wildcards: config["cov_bins"][wildcards.platform],
        plot_type=config["plot_type"],
    script:
        "../scripts/plot_precision_recall.py"


rule compare_all_tools:
    input:
        tools=lambda wildcards: expand(
            "results/{{platform}}/{{protocol}}/result_files/{method}.parquet",
            method=config["ref_tools"][wildcards.platform] + ["varlo"],
        ),
    output:
        plot=report(
            "results/{platform}/{protocol}/plots/comparisions.{plot_type}",
            category=lambda wildcards: f"{wildcards.platform} - {wildcards.protocol}",
            subcategory="All Comparisions",
            labels=lambda wildcards: {"type": "general comparisions", "method": "all"},
        ),
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/{platform}/{protocol}/compare_all_tools_{plot_type}.log",
    params:
        plot_type=config["plot_type"],
        prob_pres_threshhold=config["prob_pres_threshhold"],
    script:
        "../scripts/compare_all_tools.py"
