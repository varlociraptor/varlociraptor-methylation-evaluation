def expand_scatter_plots(platform):
    base_path = Path("results") / platform
    protocols = list(config["data"][platform].keys())
    return [str(base_path / protocol / "scatter_plot_tv.png") for protocol in protocols]


def get_protocol_sra(wildcards):
    base_path = Path("resources") / wildcards.platform / wildcards.protocol
    accession_numbers = config["data"][wildcards.platform][wildcards.protocol]
    return [
        str(base_path / SRA / "alignment_focused_dedup.bam")
        for SRA in accession_numbers
    ]
