def compute_results():
    needed_inputs = []
    # seq_platforms = ["Nanopore",  "PacBio"]
    seq_platforms = config["seq_platforms"].keys()
    for seq_platform in seq_platforms:
        # needed_inputs.append(plots(seq_platform))
        # needed_inputs.append(diff_plots(seq_platform))
        # needed_inputs.append(precision_recall(seq_platform))

        # needed_inputs.append(comparision_plots_tools(seq_platform))

        needed_inputs.append(heatmap_replicates(seq_platform))
        needed_inputs.append(correlation_table(seq_platform))
        needed_inputs.append(
            f"results/Illumina_pe/plots/heatmap_all_protocols.{config['plot_type']}"
        )

    return needed_inputs


def plots(seq_platform):
    base_path = Path("results") / seq_platform
    protocols = list(config["data"][seq_platform].keys())
    ref_methods = config["ref_tools"][seq_platform]
    return [
        str(
            base_path
            / protocol
            / (
                "plots/"
                + str(method)
                + "/plots_"
                + str(bin)
                + "."
                + config["plot_type"]
            )
        )
        for protocol in protocols
        # for method in list(ref_methods)
        for method in list(ref_methods + ["varlo"])
        for bin in list(range(0, config["cov_bins"][seq_platform])) + ["all"]
    ]


def diff_plots(seq_platform):
    base_path = Path("results") / seq_platform
    protocols = list(config["data"][seq_platform].keys())
    ref_methods = config["ref_tools"][seq_platform]
    return [
        str(
            base_path
            / protocol
            / ("plots/scatter_comp_" + method + "." + config["plot_type"])
        )
        for protocol in protocols
        for method in ref_methods
    ]


def comparision_plots_tools(seq_platform):
    base_path = Path("results") / seq_platform
    protocols = list(config["data"][seq_platform].keys())
    ref_methods = config["ref_tools"][seq_platform]
    return [
        str(base_path / protocol / ("plots/comparisions." + config["plot_type"]))
        for protocol in protocols
    ]


def heatmap_replicates(seq_platform):
    base_path = Path("results") / seq_platform
    plot_type = config["plot_type"]
    return [
        f"{base_path}/plots/{protocol}_heatmap.{plot_type}"
        for protocol in config["protocols"][seq_platform]
    ]


def correlation_table(seq_platform):
    base_path = Path("results") / seq_platform
    return [str(base_path / ("plots/correlation_table." + config["plot_type"]))]


def precision_recall(seq_platform):
    base_path = Path("results") / seq_platform
    protocols = list(config["data"][seq_platform].keys())
    return [
        f"{base_path}/{protocol}/plots/precall.{config['plot_type']}"
        for protocol in protocols
    ]


def get_protocol_sra(wildcards):
    base_path = Path("resources") / wildcards.seq_platform / wildcards.protocol
    accession_numbers = config["data"][wildcards.seq_platform][wildcards.protocol]
    return [
        str(base_path / SRA / "alignment_focused_dedup.bam")
        for SRA in accession_numbers
    ]


def get_protocol_sra_bismark(wildcards):
    base_path = Path("resources") / "ref_tools/bismark/alignment" / wildcards.protocol
    accession_numbers = config["data"]["Illumina_pe"][wildcards.protocol]
    return [
        str(base_path / SRA / (SRA + "_1_bismark_bt2_pe.bam"))
        for SRA in accession_numbers
    ]


def get_ref_methods(wildcards):
    base_path = Path("results") / "ref_tools/"
    ref_methods = config["ref_tools"][wildcards.seq_platform]
    return [
        str(base_path / method / wildcards.protocol / (method + ".bed"))
        for method in ref_methods
    ]


def get_precision_recall_csvs(wildcards):
    base_path = Path("results") / wildcards.seq_platform / wildcards.protocol / "plots"
    # methods = config["ref_tools"][wildcards.seq_platform]
    methods = config["ref_tools"][wildcards.seq_platform] + ["varlo"]
    cov_bins = [str(i) for i in range(config["cov_bins"][wildcards.seq_platform])] + [
        "all"
    ]

    return [
        str(base_path / method / f"precall_{cov_bin}.csv")
        for method in methods
        for cov_bin in cov_bins
    ]
