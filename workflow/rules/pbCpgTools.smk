# TODO MissingOutputException, even if the output is there. You have to rerun the pipeline afterwards
rule download_pb_CpG_tools:
    output:
        "resources/ref_tools/pb-CpG-tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores",
        "resources/ref_tools/pb-CpG-tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite",
    log:
        "../logs/download_pb-CpG-tools.log",
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p resources/ref_tools/pb-CpG-tools
        cd resources/ref_tools/pb-CpG-tools
        wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.1/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu.tar.gz
        tar -xzf pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu.tar.gz
        """


rule pb_CpG_tools:
    input:
        runner="resources/ref_tools/pb-CpG-tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores",
        model="resources/ref_tools/pb-CpG-tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite",
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
        {input.runner} \
        --bam {input.alignment} \
        --output-prefix {params.prefix} \
        --model {input.model} \
        --threads {threads}
        """


rule rename_pb_output:
    input:
        "results/PacBio/{protocol}/result_files/alignments_CpG.combined.bed",
    output:
        "results/PacBio/{protocol}/result_files/pb_CpG_tools.bed",
    shell:
        "mv {input} {output}"
