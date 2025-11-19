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


# This rule combines the heatmaps from  two different FDR levels into one file for easier comparison in the paper
rule heatmaps_paper:
    input:
        lambda wildcards: [
            f"results/single_sample/{wildcards.platform}/1.0/plots/{'all_samples_heatmap.parquet' if wildcards.platform == 'Illumina_pe' else 'REP_heatmap.parquet'}",
            f"results/single_sample/{wildcards.platform}/0.01/plots/{'all_samples_heatmap.parquet' if wildcards.platform == 'Illumina_pe' else 'REP_heatmap.parquet'}",
        ],
    output:
        "results/single_sample/paper/{platform}/heatmap.{plot_type}",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/heatmaps_paper/{platform}_{plot_type}.log",
    params:
        plot_type="heatmap",
        platform=lambda wildcards: wildcards.platform,
    script:
        "../scripts/paper_plots.py"


rule illumina_single_sample_paper:
    input:
        "results/single_sample/Illumina_pe/1.0/plots/bar_plot_single_samples.parquet",
        "results/single_sample/Illumina_pe/0.01/plots/bar_plot_single_samples.parquet",
    output:
        "results/single_sample/paper/Illumina_pe/illumina_single_sample.{plot_type}",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/illumina_single_sample_paper/{plot_type}.log",
    params:
        plot_type="single",
        platform="Illumina_pe",
    script:
        "../scripts/paper_plots.py"


rule bias_paper:
    input:
        "results/single_sample/Illumina_pe/0.01/plots/all_samples_bias.parquet",
        "results/single_sample/PacBio/0.01/plots/REP_bias.parquet",
        "results/single_sample/Nanopore/0.01/plots/REP_bias.parquet",
    output:
        "results/single_sample/paper/bias.{plot_type}",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/bias_paper/{plot_type}.log",
    params:
        plot_type="bias",
    script:
        "../scripts/paper_plots.py"
