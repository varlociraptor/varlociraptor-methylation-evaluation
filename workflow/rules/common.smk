def expand_scatter_plots(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    return [str(base_path / protocol / "scatter_plot_rt.png") for protocol in protocols]


def ref_tools_illumina():
    base_path = Path("results/Illumina_pe")
    protocols = list(config["data"]["Illumina_pe"].keys())
    return [
        str(base_path / protocol / "ref_tools_completed.txt") for protocol in protocols
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
