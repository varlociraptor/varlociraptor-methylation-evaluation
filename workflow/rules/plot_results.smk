rule ref_df:
    input:
        tool="results/{platform}/{protocol}/result_files/{method}.bed",
        true_meth=lambda wildcards: expand(
            "resources/bed_avg_{chrom}.bedGraph",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
    output:
        "results/{platform}/{protocol}/result_files/{method}.parquet",
    conda:
        "../envs/plot.yaml"
    params:
        cov_bin_size=lambda wildcards: config["cov_bin_size"][wildcards.platform],
        cov_bins=lambda wildcards: config["cov_bins"][wildcards.platform],
    script:
        "../scripts/df_from_calls.py"


rule plot_results_cov_specific:
    input:
        tool="results/{platform}/{protocol}/result_files/{method}.parquet",
    output:
        plot=report(
            "results/{platform}/{protocol}/plots/{method}/scatter_{cov_bin, [0-9]+}.{plot_type}",
            category=lambda wildcards: wildcards.platform,
            subcategory=lambda wildcards: f"{wildcards.method} plots",
            labels= lambda wildcards: {
                "type": "scatter plot",
                "coverage": f"{int(wildcards.cov_bin) * int(config['cov_bin_size'][wildcards.platform])} - {int(wildcards.cov_bin) * int(config['cov_bin_size'][wildcards.platform]) + int(config['cov_bin_size'][wildcards.platform])}",
                },
            ),
        tool_dist=report(
            "results/{platform}/{protocol}/plots/{method}/dist_{cov_bin, [0-9]+}.{plot_type}",
            category=lambda wildcards: wildcards.platform,
            subcategory=lambda wildcards: f"{wildcards.method} plots",
            labels= lambda wildcards: {
                "type": "dist plot",
                "coverage": f"{int(wildcards.cov_bin) * int(config['cov_bin_size'][wildcards.platform])} - {int(wildcards.cov_bin) * int(config['cov_bin_size'][wildcards.platform]) + int(config['cov_bin_size'][wildcards.platform])}",
            },
        ),
    # wildcard_constraints:
        # cov_bin=[A-Z]+
        # cov_bin="(?!all)"
    conda: 
        "../envs/plot.yaml",
    params:
        plot_type=config["plot_type"],
        prob_pres_threshhold=config["prob_pres_threshhold"],
        cov_bin=lambda wildcards: wildcards.cov_bin
    script:
        "../scripts/scatter_plot.py"

rule plot_results_all_cov:
    input:
        tool="results/{platform}/{protocol}/result_files/{method}.parquet",
    output:
        plot=report(
            "results/{platform}/{protocol}/plots/{method}/scatter_all.{plot_type}",
            # "results/{platform}/{protocol}/plots/{method}/all_scatter.{plot_type}",
            category=lambda wildcards: wildcards.platform,
            subcategory=lambda wildcards: f"{wildcards.method} plots",
            labels= lambda wildcards: {
                "type": "scatter plot",
                "coverage": "all",
                },
            ),
        tool_dist=report(
            # "results/{platform}/{protocol}/plots/{method}/dist_all.{plot_type}",
            "results/{platform}/{protocol}/plots/{method}/dist_all.{plot_type}",
            category=lambda wildcards: wildcards.platform,
            subcategory=lambda wildcards: f"{wildcards.method} plots",
            labels= lambda wildcards: {
                "type": "dist plot",
                "coverage": "all"},
        ),
    conda: 
        "../envs/plot.yaml",
    params:
        plot_type=config["plot_type"],
        prob_pres_threshhold=config["prob_pres_threshhold"],
        cov_bin = -1
    script:
        "../scripts/scatter_plot.py"

rule plot_dist_comparision:
    input:
        varlo="results/{platform}/{protocol}/result_files/varlo.parquet",
        ref_tool="results/{platform}/{protocol}/result_files/{method}.parquet"
        # ref_tool="results/ref_tools/methylDackel/TruSeq_HG002_LAB01_REP01/methylDackel.bed",
    output:
        plot=report(
            expand(
            "results/{{platform}}/{{protocol}}/plots/dist_comp_{{method}}.{plot_type}",
            plot_type=config["plot_type"],
        ),
            category=lambda wildcards: wildcards.platform,
            subcategory="All Comparisions",
            labels=lambda wildcards: {
                "type": "distance comparision",
                "method": wildcards.method,
                "coverage": "all"
            },
        ),
        scatter_plot=report(expand(
            "results/{{platform}}/{{protocol}}/plots/scatter_comp_{{method}}.{plot_type}",
            plot_type=config["plot_type"],
        ),
            category=lambda wildcards: wildcards.platform,
            subcategory="All Comparisions",
            labels=lambda wildcards: {
                "type": "scatter comparision",
                "method": wildcards.method,
                "coverage": "all"
            },
        ),
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_dist_comparision_{platform}_{protocol}_{method}.log",
    params:
        plot_type=config["plot_type"],
        cov_bins=config["cov_bins"],
        prob_pres_threshhold=config["prob_pres_threshhold"],
    script:
        "../scripts/plot_distances_together.py"


rule compute_precision_recall:
    input:
        tool="results/{platform}/{protocol}/result_files/{method}.parquet",
    output:
        precall="results/{platform}/{protocol}/plots/{method}/precall_{cov_bin}.csv",

    conda: 
        "../envs/plot.yaml",
    params:
        plot_type=config["plot_type"],
        cov_bin= lambda wildcards: wildcards.cov_bin,
        cov_bin_size=lambda wildcards: config["cov_bin_size"][wildcards.platform],
        prob_pres_threshhold=config["prob_pres_threshhold"],
    script:
        "../scripts/compute_precision_recall.py"

rule plot_precision_recall:
    input:
        get_precision_recall_csvs
        # "results/{platform}/{protocol}/plots/varlo/precall_{cov_bin}.csv",
    output:
        report(
            expand(
            "results/{{platform}}/{{protocol}}/plots/precall_{{cov_bin}}.{{plot_type}}",
            ),
            category=lambda wildcards: wildcards.platform,
            subcategory="All Comparisions",
            labels=lambda wildcards: {
                "type": "precision recall",
                "coverage": f"{int(wildcards.cov_bin) * int(config['cov_bin_size'][wildcards.platform])} - {int(wildcards.cov_bin) * int(config['cov_bin_size'][wildcards.platform]) + int(config['cov_bin_size'][wildcards.platform]) - 1}",
            },
        ),

    conda:
        "../envs/plot.yaml"
    log:
        # "logs/plot_precision_recall_{platform}_{protocol}_{plot_type}.log",
    params:
        plot_type=config["plot_type"],
    script:
        "../scripts/plot_precision_recall.py"


