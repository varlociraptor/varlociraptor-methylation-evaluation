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
        alphas=[str(config["plot_bias_alpha"])],
    script:
        "../scripts/concat_bias.py"


# Bias plot for all ARTIFACT FDR levels (one row per level)
use rule concat_plots_bias as concat_plots_bias_all_levels with:
    output:
        report(
            "results/supplementary/single_sample/bias_all_levels.pdf",
            category="supplementary",
            labels={"file": "Bias (all ARTIFACT FDR levels)"},
            caption="../report/bias.rst",
        ),
    log:
        "logs/plot_results/concat_plots_bias_all_levels.log",
    params:
        alphas=[str(a) for a in config["fdr_alpha"]],


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
