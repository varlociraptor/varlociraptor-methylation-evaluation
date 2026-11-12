rule compute_pandas_df:
    input:
        tool="results/{call_type}/{seq_platform}/called/{sample}/result_files/{method}.bed",
    output:
        "results/{call_type}/{seq_platform}/called/{sample}/result_files/{method}.parquet",
    log:
        "logs/plot_results/compute_pandas_df/{call_type}_{seq_platform}_{sample}_{method}.log",
    wildcard_constraints:
        method="(?!varlo|sample_df).*",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=34000,
    script:
        "../scripts/pandas_df_from_meth_output.py"


rule compute_varlo_df:
    input:
        # Calls already filtered by varlociraptor's FDR control at level fdr (rule control_fdr_calls)
        tool="results/{call_type}/{seq_platform}/called/{sample}/controlled_{fdr}/result_files/varlo.bed",
    output:
        "results/{call_type}/{seq_platform}/called/{sample}/result_files/varlo_{fdr}.parquet",
    log:
        "logs/plot_results/compute_varlo_df/{call_type}_{seq_platform}_{fdr}_{sample}.log",
    conda:
        "../envs/python.yaml"
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
        varlo=expand(
            "results/{{call_type}}/{{seq_platform}}/called/{{sample}}/result_files/varlo_{fdr}.parquet",
            fdr=config["fdr_alpha"],
        ),
    output:
        sample_df="results/{call_type}/{seq_platform}/result_files/sample_df_{sample}.parquet",
    log:
        "logs/plot_results/common_tool_df/{call_type}_{seq_platform}_{sample}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=60000,
    params:
        plot_type=config["plot_type"],
        filter_chrom=lambda wildcards: (
            config.get("chrom_filter", None)
            if wildcards.seq_platform == "Simulate"
            else None
        ),
    script:
        "../scripts/common_tool_df.py"


rule merge_replicates:
    input:
        lambda wildcards: expand(
            "results/{{call_type}}/{{seq_platform}}/result_files/sample_df_{rep}.parquet",
            rep=replicate_names(
                wildcards.seq_platform, wildcards.sample, wildcards.call_type
            ),
        ),
    output:
        "results/{call_type}/{seq_platform}/result_files/{sample}_replicates.parquet",
    log:
        "logs/plot_results/merge_replicates/{call_type}_{seq_platform}_{sample}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=4000,
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + ["varlo"],
    script:
        "../scripts/merge_replicates.py"


# Bin the methylation values for heatmaps and compute the discordance measurements for each tool and sample, which is used for the bias plot. This is done in one step to avoid redundant computations (binning is needed for both heatmap and bias plot).
# We offer the option to have {sample} == "all_samples", which means that the resulting df dfs will contain all samples of the respective platform. This is used for the combined heatmap and bias plot for all Illumina samples.
rule prepare_plot_df:
    input:
        "results/{call_type}/{seq_platform}/result_files/{sample}_replicates.parquet",
    output:
        df="results/{call_type}/{seq_platform}/result_files/{sample}_prepared.parquet",
        distances="results/{call_type}/{seq_platform}/result_files/{sample}_distances.parquet",
        distance_plot="results/{call_type}/{seq_platform}/result_files/{sample}_distances.html",
    log:
        "logs/plot_results/prepare_plot_df/{call_type}_{seq_platform}_{sample}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=32000,
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + [f"varlo_{fdr}" for fdr in config["fdr_alpha"]],
        sample=lambda wildcards: (
            config["samples"].get("Illumina_pe", [])
            if wildcards.sample == "all_samples"
            else [wildcards.sample]
        ),
        sample_name=lambda wildcards: wildcards.sample,
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        platform=lambda wildcards: wildcards.seq_platform,
    script:
        "../scripts/prepare_plot_df.py"
