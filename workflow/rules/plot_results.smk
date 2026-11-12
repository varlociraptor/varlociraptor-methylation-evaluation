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


rule plot_heatmaps:
    input:
        df="results/{call_type}/{seq_platform}/result_files/{sample}_prepared.parquet",
        distances="results/{call_type}/{seq_platform}/result_files/{sample}_distances.parquet",
    output:
        heatmap=report(
            "results/{call_type}/{seq_platform}/plots/{sample}_heatmap.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: seq_platform_to_name[wildcards.seq_platform],
            labels=lambda wildcards: {
                "file_type": "heatmap",
                "sample": (
                    "all samples"
                    if wildcards.sample in ["dummy", "all_samples"]
                    else wildcards.sample
                ),
            },
            caption="../report/heatmap.rst",
        ),
    # wildcard_constraints:
    #     # sample should not be 'simulated_data'
    #     sample="^(?!simulated_data$).*"
    log:
        "logs/plot_results/plot_heatmaps/{call_type}_{seq_platform}_{sample}_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=4000,
    params:
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        plot_type=lambda wildcards: wildcards.plot_type,
        callers=plot_callers,
        alpha_in_title=True,
    script:
        "../scripts/plot_heatmaps.py"


rule plot_bias:
    input:
        table="results/{call_type}/{seq_platform}/result_files/{sample}_replicates.parquet",
        # Calls passing the FDR control on the ARTIFACT event
        artifact=lambda wildcards: expand(
            "results/{call_type}/{seq_platform}/called/{rep}/result_files/varlo_artifact_{fdr}.parquet",
            call_type=wildcards.call_type,
            seq_platform=wildcards.seq_platform,
            rep=replicate_names(wildcards.seq_platform, wildcards.sample, wildcards.call_type),
            fdr=config["fdr_alpha"],
        ),
    output:
        "results/{call_type}/{seq_platform}/plots/{sample}_bias.{plot_type}",
        bias_df="results/{call_type}/{seq_platform}/plots/{sample}_bias_df_{plot_type}.parquet",
    log:
        "logs/plot_results/plot_bias/{call_type}_{seq_platform}_{sample}_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=32000,
    params:
        sample=lambda wildcards: (
            config["samples"].get("Illumina_pe", [])
            if wildcards.sample == "all_samples"
            else wildcards.sample
        ),
        platform=lambda wildcards: wildcards.seq_platform,
        plot_type=lambda wildcards: wildcards.plot_type,
        # Only the level of the combined bias plot (concat_plots_bias)
        fdrs=[config["plot_bias_alpha"]],
    script:
        "../scripts/plot_bias.py"


# Discordance vs. number of retained loci for the Varlociraptor FDR levels. Other callers get the min coverage threshold that retains the closest number of loci.
rule plot_discordance_vs_loci:
    input:
        meth_data="results/{call_type}/{seq_platform}/result_files/all_samples_replicates.parquet",
        coverages_rep1=lambda wildcards: expand(
            "results/single_sample/{{seq_platform}}/coverages/{sample}_REP01_all.regions.bed.gz",
            sample=discordance_vs_loci_samples(wildcards),
        ),
        coverages_rep2=lambda wildcards: expand(
            "results/single_sample/{{seq_platform}}/coverages/{sample}_REP02_all.regions.bed.gz",
            sample=discordance_vs_loci_samples(wildcards),
        ),
    output:
        report(
            "results/{call_type}/{seq_platform}/plots/discordance_vs_loci.{plot_type}",
            category="{call_type}",
            subcategory=lambda wildcards: seq_platform_to_name[wildcards.seq_platform],
            labels={
                "file_type": "discordance_vs_loci",
                "sample": "all samples",
            },
            caption="../report/discordance_vs_loci.rst",
        ),
        "results/{call_type}/{seq_platform}/plots/discordance_vs_loci_{plot_type}.parquet",
    log:
        "logs/plot_results/plot_discordance_vs_loci/{call_type}_{seq_platform}_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=32000,
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + [f"varlo_{fdr}" for fdr in config["fdr_alpha"]],
        samples=discordance_vs_loci_samples,
        plot_type=lambda wildcards: wildcards.plot_type,
        platform=lambda wildcards: seq_platform_to_name[wildcards.seq_platform],
    script:
        "../scripts/plot_discordance_vs_loci.py"


rule plot_runtime_comparison:
    input:
        benchmarks="benchmarks",
    output:
        tools="results/single_sample/plots/runtime_memory.{plot_type}",
    log:
        "logs/plot_results/plot_runtime_comparison/{plot_type}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_runtime_comparison.py"


rule debug_pb_cpg_tools:
    input:
        pacbio="results/single_sample/PacBio/result_files/all_samples_replicates.parquet",
        nanopore="results/single_sample/Nanopore/result_files/all_samples_replicates.parquet",
        illumina="results/single_sample/Illumina_pe/result_files/all_samples_replicates.parquet",
        pacbio_coverages_rep1=expand(
            "results/single_sample/PacBio/coverages/{sample}_REP01_all.regions.bed.gz",
            sample=config["samples"].get("PacBio", []),
        ),
        pacbio_coverages_rep2=expand(
            "results/single_sample/PacBio/coverages/{sample}_REP02_all.regions.bed.gz",
            sample=config["samples"].get("PacBio", []),
        ),
        illumina_coverages_rep1=expand(
            "results/single_sample/Illumina_pe/coverages/{sample}_REP01_all.regions.bed.gz",
            sample=config["samples"].get("Illumina_pe", []),
        ),
        illumina_coverages_rep2=expand(
            "results/single_sample/Illumina_pe/coverages/{sample}_REP02_all.regions.bed.gz",
            sample=config["samples"].get("Illumina_pe", []),
        ),
        nanopore_coverages_rep1=expand(
            "results/single_sample/Nanopore/coverages/{sample}_REP01_all.regions.bed.gz",
            sample=config["samples"].get("Nanopore", []),
        ),
        nanopore_coverages_rep2=expand(
            "results/single_sample/Nanopore/coverages/{sample}_REP02_all.regions.bed.gz",
            sample=config["samples"].get("Nanopore", []),
        ),
    output:
        report(
            "results/single_sample/plots/debug_pb_cpg_tools.{plot_type}",
            category="single_sample",
            subcategory="PacBio",
            labels={
                "file_type": "debug_pb_cpg_tools",
                "sample": "dummy",
            },
            caption="../report/debug_pb_cpg_tools.rst",
        ),
    log:
        "logs/plot_results/debug_pb_cpg_tools_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    params:
        pacbio_samples=config["samples"].get("PacBio", []),
        illumina_samples=config["samples"].get("Illumina_pe", []),
        nanopore_samples=config["samples"].get("Nanopore", []),
        # Minimum mosdepth coverage of a replicate value, applied to all callers and platforms
        min_coverage=20,
    script:
        "../scripts/debug_pb_cpg_tools.py"
