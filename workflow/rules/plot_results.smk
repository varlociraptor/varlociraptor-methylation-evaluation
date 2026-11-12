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


# We show in the paper the heatmaps only for specific FDR levels. We use this rule to compute heatmaps only for this sppecific FDR level in order to concatenate them later for the paper_plots
use rule plot_heatmaps as plot_heatmaps_paper with:
    output:
        heatmap="results/{call_type}/{seq_platform}/plots/{sample}_heatmap_paper.{plot_type}",
    log:
        "logs/plot_results/plot_heatmaps_paper/{call_type}_{seq_platform}_{sample}_{plot_type}.log",
    params:
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        plot_type=lambda wildcards: wildcards.plot_type,
        callers=lambda wildcards: plot_callers(wildcards, paper=True),
        alpha_in_title=False,


# Paper heatmaps of the single-sample long-read platforms (Modkit, pb-CpG-tools, Varlociraptor at plot_fdr_alpha)
use rule plot_heatmaps_paper as plot_heatmaps_paper_single_sample with:
    output:
        heatmap=report(
            "results/plots_paper/{call_type}/{seq_platform}/{sample}_heatmap.{plot_type}",
            category="supplementary",
            subcategory="plot_paper",
            labels=lambda wildcards: {
                "file": f"Heatmap {seq_platform_to_name[wildcards.seq_platform]}"
            },
            caption="../report/heatmap.rst",
        ),
    wildcard_constraints:
        call_type="single_sample",
        seq_platform="PacBio|Nanopore",
    log:
        "logs/plot_results/plot_heatmaps_paper_single_sample/{call_type}_{seq_platform}_{sample}_{plot_type}.log",


rule plot_bars_illumina:
    input:
        df=expand(
            "results/single_sample/Illumina_pe/result_files/{sample}_prepared.parquet",
            sample=config["samples"].get("Illumina_pe", []),
        ),
        distances=expand(
            "results/single_sample/Illumina_pe/result_files/{sample}_distances.parquet",
            sample=config["samples"].get("Illumina_pe", []),
        ),
    output:
        report(
            "results/single_sample/Illumina_pe/plots/bar_plot_single_samples.{plot_type}",
            category="supplementary",
            labels={
                "file": "Results on all Illumina samples",
            },
            caption="../report/bar_plot_illumina.rst",
        ),
    log:
        "logs/plot_results/plots_bars_illumina/single_sample_Illumina_pe_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=4000,
    params:
        sample=config["samples"].get("Illumina_pe", []),
        bin_size=lambda wildcards: config["heatmap_bin_size"],
        plot_type=lambda wildcards: wildcards.plot_type,
    script:
        "../scripts/plot_bars_illumina.py"


rule plot_bias:
    input:
        table="results/{call_type}/{seq_platform}/result_files/{sample}_replicates.parquet",
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


# Same as plot_discordance_vs_loci, but with all samples of the platform (supplementary)
use rule plot_discordance_vs_loci as plot_discordance_vs_loci_all_samples with:
    input:
        meth_data="results/{call_type}/{seq_platform}/result_files/all_samples_replicates.parquet",
        coverages_rep1=lambda wildcards: expand(
            "results/single_sample/{{seq_platform}}/coverages/{sample}_REP01_all.regions.bed.gz",
            sample=config["samples"][wildcards.seq_platform],
        ),
        coverages_rep2=lambda wildcards: expand(
            "results/single_sample/{{seq_platform}}/coverages/{sample}_REP02_all.regions.bed.gz",
            sample=config["samples"][wildcards.seq_platform],
        ),
    output:
        "results/{call_type}/{seq_platform}/plots/all_samples/discordance_vs_loci.{plot_type}",
        "results/{call_type}/{seq_platform}/plots/all_samples/discordance_vs_loci_{plot_type}.parquet",
    log:
        "logs/plot_results/plot_discordance_vs_loci_all_samples/{call_type}_{seq_platform}_{plot_type}.log",
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get(
            wildcards.seq_platform, []
        )
        + [f"varlo_{fdr}" for fdr in config["fdr_alpha"]],
        samples=lambda wildcards: config["samples"][wildcards.seq_platform],
        plot_type=lambda wildcards: wildcards.plot_type,
        platform=lambda wildcards: seq_platform_to_name[wildcards.seq_platform],


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


rule concat_plots_multi_sample:
    input:
        plots=[
            "results/multi_sample/pb_methylSeq/plots/dummy_heatmap_paper.pdf",
            "results/multi_sample/np_pb/plots/dummy_heatmap_paper.pdf",
            "results/multi_sample/np_methylSeq/plots/dummy_heatmap_paper.pdf",
        ],
    output:
        report(
            "results/plots_paper/multi_sample/combined_heatmap.pdf",
            category="supplementary",
            subcategory="plot_paper",
            labels={"file": "Multi-sample heatmaps"},
            caption="../report/combined_heatmap_multi_sample.rst",
        ),
    log:
        "logs/plot_results/concat_plots_multi_sample.log",
    conda:
        "../envs/fitz.yaml"
    resources:
        mem_mb=4000,
    script:
        "../scripts/concat_plots.py"


rule concat_plots_bias:
    input:
        illumina="results/single_sample/Illumina_pe/plots/all_samples_bias_df_pdf.parquet",
        pacbio="results/single_sample/PacBio/plots/dummy_bias_df_pdf.parquet",
        nanopore="results/single_sample/Nanopore/plots/dummy_bias_df_pdf.parquet",
    output:
        report(
            "results/plots_paper/single_sample/bias.pdf",
            category="supplementary",
            subcategory="plot_paper",
            labels={"file": "bias"},
            caption="../report/bias.rst",
        ),
    log:
        "logs/plot_results/concat_plots_bias.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=4000,
    params:
        alpha=str(config["plot_bias_alpha"]),
    script:
        "../scripts/concat_bias.py"


rule concat_plots_coverage:
    input:
        illumina="results/single_sample/Illumina_pe/coverages/all_samples_coverage_plot_pdf_all.parquet",
        pacbio="results/single_sample/PacBio/coverages/dummy_coverage_plot_pdf_all.parquet",
        nanopore="results/single_sample/Nanopore/coverages/dummy_coverage_plot_pdf_all.parquet",
    output:
        report(
            "results/plots_paper/single_sample/coverage.pdf",
            category="supplementary",
            subcategory="plot_paper",
            labels={"file": "Coverage retention"},
            caption="../report/coverage_retained.rst",
        ),
    log:
        "logs/plot_results/concat_plots_coverage.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=4000,
    params:
        callers=[
            tool for tools in config["ref_tools"].values() for tool in tools
        ]
        + [f"varlo_{fdr}" for fdr in config["plot_fdr_alpha"].get("coverage", config["fdr_alpha"])],
    script:
        "../scripts/concat_coverage.py"


rule concat_plots_discordance_vs_loci:
    input:
        illumina="results/{call_type}/Illumina_pe/plots/discordance_vs_loci_{plot_type}.parquet",
        pacbio="results/{call_type}/PacBio/plots/discordance_vs_loci_{plot_type}.parquet",
        nanopore="results/{call_type}/Nanopore/plots/discordance_vs_loci_{plot_type}.parquet",
    output:
        report(
            "results/plots_paper/{call_type}/discordance_vs_loci.{plot_type}",
            category="supplementary",
            subcategory="plot_paper",
            labels={"file": "Discordance vs. loci"},
            caption="../report/discordance_vs_loci.rst",
        ),
    log:
        "logs/plot_results/concat_plots_discordance_vs_loci/{call_type}_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=4000,
    script:
        "../scripts/concat_discordance_vs_loci.py"


# Discordance vs. loci with all samples of each platform (supplementary)
use rule concat_plots_discordance_vs_loci as concat_plots_discordance_vs_loci_all_samples with:
    input:
        illumina=lambda wildcards: discordance_vs_loci_all_parquet(
            wildcards.call_type, "Illumina_pe", wildcards.plot_type
        ),
        pacbio=lambda wildcards: discordance_vs_loci_all_parquet(
            wildcards.call_type, "PacBio", wildcards.plot_type
        ),
        nanopore=lambda wildcards: discordance_vs_loci_all_parquet(
            wildcards.call_type, "Nanopore", wildcards.plot_type
        ),
    output:
        report(
            "results/supplementary/{call_type}/discordance_vs_loci_all_samples.{plot_type}",
            category="supplementary",
            labels={"file": "Discordance vs. loci (all samples)"},
            caption="../report/discordance_vs_loci.rst",
        ),
    log:
        "logs/plot_results/concat_plots_discordance_vs_loci_all_samples/{call_type}_{plot_type}.log",


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


# Static example of the formula tree that Varlociraptor traverses to evaluate an event
rule event_tree:
    input:
        workflow.source_path("../resources/event_tree.pdf"),
    output:
        report(
            "results/supplementary/event_tree.pdf",
            category="supplementary",
            labels={
                "file": f"Exemplary event_tree",
            },
            caption="../report/event_tree.rst",
        ),
    log:
        "logs/plot_results/event_tree.log",
    shell:
        "cp {input} {output} 2> {log}"
