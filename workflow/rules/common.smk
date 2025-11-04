from pathlib import Path
from typing import List, Union


def compute_results() -> List[List[str]]:
    """
    Collect all input file paths required for the workflow.
    """
    inputs: List[List[str]] = []

    # Heatmaps per sequencing platform
    for platform in config["seq_platforms"].keys():
        inputs.append(heatmap_replicates(platform))

    # Single-sample heatmaps across all FDR thresholds
    if "Illumina_pe" in config["seq_platforms"]:
        inputs.append(
            [
                f"results/single_sample/Illumina_pe/{fdr}/plots/heatmap_all_samples.{config['plot_type']}"
                for fdr in config["fdr_alpha"]
            ]
        )

    # Multi-sample common heatmaps
    inputs.append(heatmap_replicates_common())

    return inputs


def heatmap_replicates(seq_platform: str) -> List[str]:
    """
    Return file paths for replicate heatmaps for a given sequencing platform.
    """
    base_path = Path("results/single_sample") / seq_platform
    plot_type = config["plot_type"]

    return [
        f"{base_path}/{fdr}/plots/{sample}_heatmap.{plot_type}"
        for sample in config["samples"][seq_platform]
        for fdr in config["fdr_alpha"]
    ]


def heatmap_replicates_common() -> List[str]:
    """
    Return file paths for multi-sample common heatmaps across comparisons.
    """
    base_path = Path("results/multi_sample")
    plot_type = config["plot_type"]

    comparisons = ["np_pb", "pb_trueOX", "np_trueOX"]

    return [
        f"{base_path}/{comp}/{fdr}/plots/{sample}_heatmap.{plot_type}"
        for comp in comparisons
        for sample in config["samples"].get("multi_sample", [])
        for fdr in config["fdr_alpha"]
    ]


def get_sample_sra(wildcards) -> List[str]:
    """
    Return BAM file paths for a given platform and sample, based on the config.
    """
    base_path = Path("resources") / wildcards.seq_platform / wildcards.sample

    if wildcards.seq_platform not in config["data"]:
        return []
    if wildcards.sample not in config["data"].get(wildcards.seq_platform, {}):
        return []

    accession_numbers = (
        config["data"].get(wildcards.seq_platform, {}).get(wildcards.sample, [])
    )
    return [
        str(base_path / sra / "alignment_focused_dedup.bam")
        for sra in accession_numbers
    ]


def get_sample_sra_bismark(wildcards) -> List[str]:
    """
    Return Bismark alignment BAM file paths for a given sample.
    """
    base_path = Path("resources/ref_tools/bismark/alignment") / wildcards.sample
    accession_numbers = config["data"]["Illumina_pe"][wildcards.sample]

    return [
        f"resources/ref_tools/bismark/bams/{wildcards.sample}_pe_{sra}_unsorted.bam"
        for sra in accession_numbers
    ]
