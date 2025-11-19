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


def compute_all_samples(wildcards):
    tools = [
        f"results/{wildcards.call_type}/{wildcards.seq_platform}/called/{sample}/result_files/{method}.parquet"
        for sample in config["data"].get(wildcards.seq_platform, [])
        for method in config["ref_tools"].get(wildcards.seq_platform, [])
    ]
    print([sample for sample in config["data"].get(wildcards.seq_platform, [])])
    varlo = [
        f"results/{wildcards.call_type}/{wildcards.seq_platform}/{fdr}/{sample}/result_files/varlo.parquet"
        for sample in config["data"].get(wildcards.seq_platform, [])
        for fdr in config["fdr_alpha"]
    ]
    return tools + varlo


# Computes one common df out of all single methylation caller dfs
rule common_tool_df:
    input:
        tools=lambda wildcards: expand(
            "results/{{call_type}}/{{seq_platform}}/called/{{sample}}/result_files/{method}.parquet",
            method=config["ref_tools"].get(wildcards.seq_platform, []),
        ),
        varlo=expand(
            "results/{{call_type}}/{{seq_platform}}/{fdr}/{{sample}}/result_files/varlo.parquet",
            fdr=config["fdr_alpha"],
        ),
    output:
        sample_df="results/{call_type}/{seq_platform}/{sample}/result_files/sample_df.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/common_tool_df/{call_type}_{seq_platform}_{sample}.log",
    params:
        plot_type=config["plot_type"],
    resources:
        mem_mb=64000,
    script:
        "../scripts/common_tool_df.py"


rule merge_replicates:
    input:
        samples=lambda wildcards: expand(
            "results/{call_type}/{seq_platform}/{sample}/result_files/sample_df.parquet",
            call_type=wildcards.call_type,
            seq_platform=wildcards.seq_platform,
            sample=config["data"].get(wildcards.seq_platform, wildcards.seq_platform),
        ),
    output:
        table="results/{call_type}/{seq_platform}/plots/replicates.hd5",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/merge_replicates/{call_type}_{seq_platform}.log",
    # wildcard_constraints:
    #     call_type="(?!multi_sample).*",
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
        "../scripts/merge_replicates.py"


rule prepare_plot_df:
    input:
        "results/{call_type}/{seq_platform}/plots/replicates.hd5",
    output:
        df="results/{call_type}/{seq_platform}/plots/{sample}_prepared.parquet",
        mapes="results/{call_type}/{seq_platform}/plots/{sample}_mapes.parquet",
        # biases="results/{call_type}/{seq_platform}/plots/{sample}_biases.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/prepare_plot_df/{call_type}_{seq_platform}_{sample}.log",
    resources:
        mem_mb=64000,
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + [f"varlo_{fdr}" for fdr in config["fdr_alpha"]],
        sample=lambda wildcards: (
            config["samples"].get("Illumina_pe", [])
            if wildcards.sample == "all_samples"
            else wildcards.sample
        ),
        sample_name=lambda wildcards: wildcards.sample,
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        platform=lambda wildcards: wildcards.seq_platform,
    script:
        "../scripts/prepare_plot_df.py"


rule plot_heatmaps:
    input:
        df="results/{call_type}/{seq_platform}/plots/{sample}_prepared.parquet",
        mapes="results/{call_type}/{seq_platform}/plots/{sample}_mapes.parquet",
    output:
        heatmap=report(
            "results/{call_type}/{seq_platform}/plots/{sample}_heatmap.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: f"{wildcards.seq_platform}",
            labels={
                "file": "heatmap",
                "sample": "{sample}",
            },
            caption="../report/heatmap.rst",
        ),
    conda:
        "../envs/plot.yaml"
    resources:
        mem_mb=64000,
    log:
        "logs/plot_results/plot_heatmaps/{call_type}_{seq_platform}_{sample}_{plot_type}.log",
    params:
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        plot_type=lambda wildcards: wildcards.plot_type,
    script:
        "../scripts/plot_heatmaps.py"


# Compute common heatmap over all Illumina samples
rule plots_bars_illumina:
    input:
        df=expand(
            "results/single_sample/Illumina_pe/plots/{sample}_prepared.parquet",
            sample=config["samples"].get("Illumina_pe", []),
        ),
        mapes=expand(
            "results/single_sample/Illumina_pe/plots/{sample}_mapes.parquet",
            sample=config["samples"].get("Illumina_pe", []),
        ),
    output:
        report(
            "results/single_sample/Illumina_pe/plots/bar_plot_single_samples.{plot_type}",
            category="single_sample",
            subcategory="Illumina_pe",
            labels={
                "file": "bar_plot_single_samples",
                "sample": "all samples",
            },
            caption="../report/bar_plot_illumina.rst",
        ),
    conda:
        "../envs/plot.yaml"
    resources:
        mem_mb=64000,
    log:
        "logs/plot_results/plots_bars_illumina/single_sample_Illumina_pe_{plot_type}.log",
    params:
        sample=config["samples"].get("Illumina_pe", []),
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        plot_type=lambda wildcards: wildcards.plot_type,
    script:
        "../scripts/plot_bars_illumina.py"


rule plot_bias:
    input:
        "results/{call_type}/{seq_platform}/{fdr}/plots/replicates.hd5",
    output:
        report(
            "results/{call_type}/{seq_platform}/plots/{sample}_bias_{fdr}.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: f"{wildcards.seq_platform}",
            labels={
                "file": "bias",
                "sample": "{sample}",
                "fdr": "{fdr}",
            },
            caption="../report/bias.rst",
        ),
    conda:
        "../envs/plot.yaml"
    resources:
        mem_mb=64000,
    log:
        "logs/plot_results/plot_bias/{call_type}_{seq_platform}_{fdr}_{sample}_{plot_type}.log",
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + ["varlo"],
        sample=lambda wildcards: (
            config["samples"].get("Illumina_pe", [])
            if wildcards.sample == "all_samples"
            else wildcards.sample
        ),
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        fdr=lambda wildcards: wildcards.fdr,
        platform=lambda wildcards: wildcards.seq_platform,
        plot_type=lambda wildcards: wildcards.plot_type,
    script:
        "../scripts/plot_bias.py"


rule plot_runtime_comparison:
    input:
        benchmarks="benchmarks",
    output:
        tools="results/single_sample/paper/runtime_memory.{plot_type}",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/plot_runtime_comparison/{plot_type}.log",
    script:
        "../scripts/plot_runtime_comparison.py"
