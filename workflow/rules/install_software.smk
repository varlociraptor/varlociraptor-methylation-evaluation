rule install_modkit:
    output:
        "modkit_installed.flag",
    conda:
        "../envs/modkit.yaml"
    log:
        "../logs/install_modkit.log",
    shell:
        """
        export PATH=$PATH:~/.cargo/bin
        export PATH=$PATH:/homes/aprinz/.cargo/bin
        cargo install --git https://github.com/nanoporetech/modkit.git
        touch {output}
        """


rule download_pb_CpG_tools:
    output:
        directory("../../pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu"),
    log:
        "../logs/download_pb-CpG-tools.log",
    params:
        base_dir=config["base_dir"],
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        cd {params.base_dir}
        wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.1/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu.tar.gz
        tar -xzf pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu.tar.gz
        """
