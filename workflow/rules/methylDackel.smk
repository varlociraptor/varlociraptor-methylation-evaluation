# Vielleicht vorher Regel mit mkdir {params.pipeline_path}/results/ref_tools/dackel/{wildcards.protocol}/ einfuegen, konnte ich jetzt noch nicht ueberpruefen, gibt aber sonst glaub ich nen Fehler
rule methylDackel:
    input:
        genome="resources/genome.fasta",
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        "results/ref_tools/methylDackel/{protocol}/alignments_CpG.bedGraph",
    conda:
        "../envs/methylDackel.yaml"
    log:
        "logs/methylDackel_{protocol}.log",
    params:
        pipeline_path=config["pipeline_path"],
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ".bedGraph"
        ),
    shell:
        """
        MethylDackel extract {input.genome} {input.alignment} -o {params.pipeline_path}results/ref_tools/methylDackel/{wildcards.protocol}/alignments --mergeContext
        """


rule rename_dackel_output:
    input:
        "results/ref_tools/methylDackel/{protocol}/alignments_CpG.bedGraph",
    output:
        "results/ref_tools/methylDackel/{protocol}/methylDackel.bed",
    shell:
        "mv {input} {output}"
