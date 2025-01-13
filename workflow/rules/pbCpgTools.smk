rule download_pb_CpG_tools:
    output:
        directory("../tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu"),
    log:
        "../logs/download_pb-CpG-tools.log",
    params:
        base_dir=config["base_dir"],
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p {params.base_dir}/tools
        cd {params.base_dir}/tools
        wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.1/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu.tar.gz
        tar -xzf pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu.tar.gz
        """


rule pb_CpG_tools:
    input:
        pb_tools_dir="../tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu",
        alignment="resources/PacBio/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/PacBio/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        genome="resources/genome.fasta",
    output:
        "results/PacBio/{protocol}/result_files/alignments_CpG.combined.bed",
    log:
        "logs/pb_CpG_tools_{protocol}.log",
    params:
        base_dir=config["base_dir"],
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ""
        ),
    threads: 8
    shell:
        """
        {params.base_dir}pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
        --bam {input.alignment} \
        --output-prefix {params.prefix} \
        --model {params.base_dir}pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
        --threads {threads}
        """


rule rename_pb_output:
    input:
        "results/PacBio/{protocol}/result_files/alignments_CpG.combined.bed",
    output:
        "results/PacBio/{protocol}/result_files/pb_CpG_tools.bed",
    shell:
        "mv {input} {output}"
