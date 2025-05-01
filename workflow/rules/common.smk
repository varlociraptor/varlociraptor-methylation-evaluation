def compute_results():
    needed_inputs = []
    for platform in config["platforms"].keys():
        needed_inputs.append(plots(platform))
        needed_inputs.append(diff_plots(platform))
        needed_inputs.append(precision_recall(platform))
        needed_inputs.append(comparision_plots_tools(platform))
        needed_inputs.append(comparision_plots_samples(platform))

    return needed_inputs


def plots(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    ref_methods = config["ref_tools"][platform]
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
        for method in list(ref_methods + ["varlo"])
        for bin in list(range(0, config["cov_bins"][platform])) + ["all"]
    ]


def diff_plots(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    ref_methods = config["ref_tools"][platform]
    return [
        str(
            base_path
            / protocol
            / ("plots/scatter_comp_" + method + "." + config["plot_type"])
        )
        for protocol in protocols
        for method in ref_methods
    ]


def comparision_plots_tools(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    ref_methods = config["ref_tools"][platform]
    return [
        str(base_path / protocol / ("plots/comparisions." + config["plot_type"]))
        for protocol in protocols
    ]


def comparision_plots_samples(platform):
    base_path = Path("results") / platform
    ref_methods = config["ref_tools"][platform]
    return [str(base_path / ("plots/comparisions." + config["plot_type"]))]


def precision_recall(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    return [
        f"{base_path}/{protocol}/plots/precall.{config['plot_type']}"
        for protocol in protocols
    ]


def get_protocol_sra(wildcards):
    base_path = Path("resources") / wildcards.platform / wildcards.protocol
    accession_numbers = config["data"][wildcards.platform][wildcards.protocol]
    return [
        str(base_path / SRA / "alignment_focused_dedup.bam")
        for SRA in accession_numbers
    ]


def get_protocol_sra_bismark(wildcards):
    base_path = Path("resources") / "ref_tools/bismark/alignment" / wildcards.protocol
    accession_numbers = config["data"]["Illumina_pe"][wildcards.protocol]
    return [
        str(base_path / SRA / (SRA + "_1_trimmed_bismark_bt2_pe.bam"))
        for SRA in accession_numbers
    ]


def get_ref_methods(wildcards):
    base_path = Path("results") / "ref_tools/"
    ref_methods = config["ref_tools"][wildcards.platform]
    return [
        str(base_path / method / wildcards.protocol / (method + ".bed"))
        for method in ref_methods
    ]


def get_precision_recall_csvs(wildcards):
    base_path = Path("results") / wildcards.platform / wildcards.protocol / "plots"
    methods = config["ref_tools"][wildcards.platform] + ["varlo"]
    cov_bins = [str(i) for i in range(config["cov_bins"][wildcards.platform])] + ["all"]

    return [
        str(base_path / method / f"precall_{cov_bin}.csv")
        for method in methods
        for cov_bin in cov_bins
    ]
