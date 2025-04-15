# TODO MissingOutputException, even if the output is there. You have to rerun the pipeline afterwards
rule pb_CpG_download:
    output:
        "resources/ref_tools/pb-CpG-tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores",
        "resources/ref_tools/pb-CpG-tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite",
    log:
        "../logs/pb_CpG/download.log",
    conda:
        "../envs/shell_cmds.yaml"
    shell:
        """
        mkdir -p resources/ref_tools/pb-CpG-tools 2> {log}
        cd resources/ref_tools/pb-CpG-tools 2> {log}
        wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.1/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu.tar.gz 2> {log}
        tar -xzf pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu.tar.gz 2> {log}
        """


rule pb_CpG_compute_methylation:
    input:
        runner="resources/ref_tools/pb-CpG-tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores",
        model="resources/ref_tools/pb-CpG-tools/pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite",
        alignment="resources/PacBio/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/PacBio/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        genome="resources/genome.fasta",
    output:
        "results/PacBio/{protocol}/result_files/alignments_CpG.combined.bed",
    log:
        "logs/pb_CpG/{protocol}/compute_methylation.log",
    params:
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ""
        ),
    threads: 8
    shell:
        "{input.runner} --bam {input.alignment} --output-prefix {params.prefix} --model {input.model} --threads {threads} 2> {log}"


rule pb_CpG_rename_output:
    input:
        "results/PacBio/{protocol}/result_files/alignments_CpG.combined.bed",
    output:
        "results/PacBio/{protocol}/result_files/pb_CpG_tools.bed",
    log:
        "logs/pb_CpG/{protocol}/rename_output.log",
    shell:
        "mv {input} {output} 2> {log}"
