def compute_results():
    needed_inputs = []
    for platform in config["platforms"]:
        needed_inputs.append(scatter_plots(platform))
        needed_inputs.append(scatter_plots_ref(platform))
    return needed_inputs


def scatter_plots(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    return [
        str(base_path / protocol / ("varlo_" + str(bin) + ".png"))
        for protocol in protocols
        for bin in range(0, config["cov_bins"])
    ]


def scatter_plots_ref(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    ref_methods = config["ref_tools"][platform]
    return [
        str(base_path / protocol / (method + ".png"))
        for protocol in protocols
        for method in ref_methods
    ]


def get_protocol_sra(wildcards):
    base_path = Path("resources") / wildcards.platform / wildcards.protocol
    accession_numbers = config["data"][wildcards.platform][wildcards.protocol]
    return [
        str(base_path / SRA / "alignment_focused_dedup.bam")
        for SRA in accession_numbers
    ]
    # + [str(base_path / SRA / "ref_tools_completed.txt") for SRA in accession_numbers]


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
