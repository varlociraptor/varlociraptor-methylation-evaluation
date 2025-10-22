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
    inputs.append(
        [
            f"results/single_sample/Illumina_pe/{fdr}/plots/heatmap_all_protocols.{config['plot_type']}"
            for fdr in config["fdr_alpha"]
        ]
    )

    # Multi-sample common heatmaps
    inputs.append(heatmap_replicates_common())

    # Runtime comparison plots
    inputs.append([f"results/runtime_comparison_tools.{config['plot_type']}"])
    inputs.append([f"results/runtime_comparison_varlo.{config['plot_type']}"])

    return inputs


def heatmap_replicates(seq_platform: str) -> List[str]:
    """
    Return file paths for replicate heatmaps for a given sequencing platform.
    """
    base_path = Path("results/single_sample") / seq_platform
    plot_type = config["plot_type"]

    return [
        f"{base_path}/{fdr}/plots/{protocol}_heatmap.{plot_type}"
        for protocol in config["protocols"][seq_platform]
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
        f"{base_path}/{comp}/{fdr}/plots/{protocol}_heatmap.{plot_type}"
        for comp in comparisons
        for protocol in config["protocols"]["multi_sample"]
        for fdr in config["fdr_alpha"]
    ]


def get_protocol_sra(wildcards) -> List[str]:
    """
    Return BAM file paths for a given platform and protocol, based on the config.
    """
    base_path = Path("resources") / wildcards.seq_platform / wildcards.protocol

    if wildcards.seq_platform not in config["data"]:
        return []
    if wildcards.protocol not in config["data"][wildcards.seq_platform]:
        return []

    accession_numbers = config["data"][wildcards.seq_platform][wildcards.protocol]
    return [
        str(base_path / sra / "alignment_focused_dedup.bam")
        for sra in accession_numbers
    ]


def get_protocol_sra_bismark(wildcards) -> List[str]:
    """
    Return Bismark alignment BAM file paths for a given protocol.
    """
    base_path = Path("resources/ref_tools/bismark/alignment") / wildcards.protocol
    accession_numbers = config["data"]["Illumina_pe"][wildcards.protocol]

    return [
        str(base_path / sra / f"{sra}_1_bismark_bt2_pe.bam")
        for sra in accession_numbers
    ]
