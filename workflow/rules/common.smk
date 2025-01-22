def compute_results():
    needed_inputs = []
    for platform in config["platforms"].keys():
        needed_inputs.append(scatter_plots(platform))
        needed_inputs.append(diff_plots(platform))
        needed_inputs.append(precision_recall(platform))
    return needed_inputs


def scatter_plots(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    ref_methods = config["ref_tools"][platform]
    return [
        # str(base_path / protocol / ("plots/" + method + str(bin) + snakemake.params["plot_type"]))
        str(
            base_path
            / protocol
            / (
                "plots/"
                + str(method)
                + "/"
                + type
                + "_"
                + str(bin)
                + "."
                + config["plot_type"]
            )
        )
        for protocol in protocols
        for type in ["scatter", "dist"]
        for method in list(ref_methods + ["varlo"])
        for bin in list(range(0, config["cov_bins"][platform])) + ["all"]
        # for bin in list(range(0, config["cov_bins"][platform]))
    ]


def diff_plots(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    ref_methods = config["ref_tools"][platform]
    return [
        str(
            base_path
            / protocol
            / ("plots/" + plot_type + "_comp_" + method + "." + config["plot_type"])
        )
        for plot_type in ["scatter", "dist"]
        for protocol in protocols
        for method in ref_methods
    ]


def precision_recall(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    return [
        f"{base_path}/{protocol}/plots/precall_{bin}.{config['plot_type']}"
        for protocol in protocols
        for bin in list(range(0, config["cov_bins"][platform])) + ["all"]
        # for bin in range(0, config["cov_bins"][platform])
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
        str(base_path / SRA / (SRA + "_1_trimmed_bismark_bt2_pe.deduplicated.bam"))
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
    methods = config["ref_tools"][wildcards.platform]
    methods.append("varlo")
    return [str(base_path / method / "precall_{cov_bin}.csv") for method in methods]
