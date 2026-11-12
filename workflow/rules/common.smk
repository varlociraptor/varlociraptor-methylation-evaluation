from pathlib import Path
from typing import List, Union


def compute_results() -> List[List[str]]:
    """
    Collect all input file paths required for the workflow.
    """
    inputs: List[List[str]] = []

    # Heatmaps per sequencing platform
    for platform in config["seq_platforms"].keys():
        if platform != "Simulate" and platform != "multi_sample":
            inputs.append(heatmap_replicates(platform))
            inputs.append(
                f"results/single_sample/{platform}/plots/discordance_vs_loci.{config['plot_type']}"
            )

    # Single-sample heatmaps across all FDR thresholds
    if "Illumina_pe" in config["seq_platforms"]:
        inputs.append(
            f"results/single_sample/Illumina_pe/plots/all_samples_heatmap.{config['plot_type']}"
        )
        inputs.append(
            f"results/single_sample/Illumina_pe/plots/bar_plot_single_samples.{config['plot_type']}"
        )
    if "Simulate" in config["seq_platforms"]:

        chromosome = (
            config["seq_platforms"]["Simulate"]
            if config["seq_platforms"]["Simulate"] != "genome"
            else "21"
        )
        inputs.append(
            f"results/single_sample/Simulate/plots/simulated_data_{chromosome}.{config['plot_type']}"
        )
        inputs.append(
            f"results/single_sample/Simulate/plots/fdr_calibration_{chromosome}.{config['plot_type']}"
        )
        inputs.append(
            f"results/single_sample/Simulate/plots/pit_calibration_{chromosome}.{config['plot_type']}"
        )
    # Multi-sample common heatmaps
    if "multi_sample" in config["seq_platforms"]:
        inputs.append(heatmap_replicates_common())
    if "PacBio" in config["seq_platforms"]:
        inputs.append(
            [
                f"results/single_sample/plots/debug_pb_cpg_tools.{config['plot_type']}",
            ]
        )

    # Supplementary: example of an event's formula tree
    inputs.append("results/supplementary/event_tree.pdf")

    # Combined paper figures (concat_* rules)
    if config["paper_plots"]:
        inputs.append(paper_plots())

    return inputs


def heatmap_replicates(seq_platform: str) -> List[str]:
    """
    Return file paths for replicate heatmaps for a given sequencing platform.
    """
    base_path = Path("results/single_sample") / seq_platform
    plot_type = config["plot_type"]

    return [
        f"{base_path}/plots/{sample}_heatmap.{plot_type}"
        for sample in config["samples"][seq_platform]
    ]


def heatmap_replicates_common() -> List[str]:
    """
    Return file paths for multi-sample common heatmaps across comparisons.
    """
    base_path = Path("results/multi_sample")
    plot_type = config["plot_type"]

    comparisons = ["np_pb", "pb_methylSeq", "np_methylSeq"]
    return [
        f"{base_path}/{comp}/plots/{sample}_heatmap.{plot_type}"
        for comp in comparisons
        for sample in config["samples"].get("multi_sample", [])
    ]


def plot_callers(wildcards, paper: bool = False) -> List[str]:
    """Heatmap callers: reference tools, then the Varlociraptor levels (paper: only those of plot_fdr_alpha)."""
    fdrs = config["fdr_alpha"]
    if paper:
        key = "multi_sample" if wildcards.call_type == "multi_sample" else wildcards.seq_platform
        fdrs = config["plot_fdr_alpha"].get(key, fdrs)
    return config["ref_tools"].get(wildcards.seq_platform, []) + [f"varlo_{fdr}" for fdr in fdrs]

def paper_plots() -> List[str]:
    """
    Return file paths of all combined paper figures (concat_* rules) under results/plots_paper/.
    """
    plot_type = config["plot_type"]
    platforms = config["seq_platforms"]
    paths = []
    if all(platform in platforms for platform in ["Illumina_pe", "PacBio", "Nanopore"]):
        paths.append("results/plots_paper/single_sample/bias.pdf")
        paths.append(
            f"results/plots_paper/single_sample/discordance_vs_loci.{plot_type}"
        )
        paths.append("results/supplementary/single_sample/bias_all_levels.pdf")
        paths.append(
            f"results/supplementary/single_sample/discordance_vs_loci_all_samples.{plot_type}"
        )
    for platform in ["PacBio", "Nanopore"]:
        if platform in platforms:
            paths.append(
                f"results/plots_paper/single_sample/{platform}/dummy_heatmap.{plot_type}"
            )
    if "multi_sample" in platforms:
        paths.append("results/plots_paper/multi_sample/combined_heatmap.pdf")
    return paths


def discordance_vs_loci_samples(wildcards) -> List[str]:
    return config.get("discordance_vs_loci_samples", {})[wildcards.seq_platform]


def discordance_vs_loci_all_parquet(
    call_type: str, platform: str, plot_type: str
) -> str:
    """
    Parquet of discordance vs. loci with all samples of the platform. Reuses the one of the selected samples
    if the selection already contains all samples (avoids recomputing).
    """
    selected = config.get("discordance_vs_loci_samples", {}).get(platform, [])
    if sorted(selected) == sorted(config["samples"].get(platform, [])):
        return f"results/{call_type}/{platform}/plots/discordance_vs_loci_{plot_type}.parquet"
    return f"results/{call_type}/{platform}/plots/all_samples/discordance_vs_loci_{plot_type}.parquet"


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
    accession_numbers = config["data"][wildcards.platform][wildcards.sample]

    return [
        f"resources/ref_tools/bismark/{wildcards.platform}/bams/{wildcards.sample}_pe_{sra}_unsorted.bam"
        for sra in accession_numbers
    ]


seq_platform_to_name = {
    "np_methylSeq": "Nanopore & MethylSeq",
    "pb_methylSeq": "PacBio & MethylSeq",
    "np_pb": "Nanopore & PacBio",
    "Illumina_pe": "Illumina",
    "Nanopore": "Nanopore",
    "PacBio": "PacBio",
}


def varlociraptor_binary():
    if config.get("varlociraptor-branch") is None:
        return []
    return "resources/tools/meth_calling/varlociraptor/target/release/varlociraptor"


OMIT_BIAS_FLAGS = [
    "--omit-strand-bias",
    "--omit-read-orientation-bias",
    "--omit-read-position-bias",
    "--omit-softclip-bias",
    "--omit-homopolymer-artifact-detection",
    "--omit-alt-locus-bias",
]


def omit_bias_flags() -> str:
    """Flags of `varlociraptor call variants` to disable the bias model (config omit_bias)."""
    return " ".join(OMIT_BIAS_FLAGS) if config.get("omit_bias", False) else ""


def varlociraptor_cmd(wildcards, input):
    if config.get("varlociraptor-branch") is None:
        return "varlociraptor"
    return f"chmod +x {input.varlo}; {input.varlo}"


def replicate_names(seq_platform, sample, call_type):
    """
    Return the replicates of a sample.
    "dummy" and "all_samples" return all replicates of the respective platform.
    """
    replicates = list(
        config["data"].get(seq_platform, config["data"].get(call_type, []))
    )
    if sample in ["all_samples", "dummy"]:
        return replicates
    return [rep for rep in replicates if rep.startswith(f"{sample}_REP")]
