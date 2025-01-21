# Vielleicht vorher Regel mit mkdir {params.pipeline_path}/results/ref_tools/dackel/{wildcards.protocol}/ einfuegen, konnte ich jetzt noch nicht ueberpruefen, gibt aber sonst glaub ich nen Fehler
rule methylDackel:
    input:
        genome="resources/genome.fasta",
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        "results/Illumina_pe/{protocol}/result_files/alignments_CpG.bedGraph",
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
        echo test
        OUTDIR=$(dirname {output})/alignments
        echo $OUTDIR
        MethylDackel extract {input.genome} {input.alignment} -o $OUTDIR --mergeContext
        """


# Vielleicht vorher Regel mit mkdir {params.pipeline_path}/results/ref_tools/dackel/{wildcards.protocol}/ einfuegen, konnte ich jetzt noch nicht ueberpruefen, gibt aber sonst glaub ich nen Fehler
rule methylDackel_debug:
    input:
        genome="resources/genome.fasta",
        # alignment="resources/Illumina_pe/TruSeq_HG002_LAB01_REP01/alignment_focused_downsampled_dedup_renamed.bam",
        # alignment_index="resources/Illumina_pe/TruSeq_HG002_LAB01_REP01/alignment_focused_downsampled_dedup_renamed.bam.bai",
        alignment="resources/Illumina_pe/TruSeq_HG002_LAB01_REP01/candidate_specific/alignment_valid_20-of-20.bam",
        alignment_index="resources/Illumina_pe/TruSeq_HG002_LAB01_REP01/candidate_specific/alignment_valid_20-of-20.bam.bai",
    output:
        "resources/Illumina_pe/TruSeq_HG002_LAB01_REP01/candidate_specific/alignments_CpG.bedGraph",
    conda:
        "../envs/methylDackel.yaml"
    log:
        "logs/methylDackel.log",
    params:
        pipeline_path=config["pipeline_path"],
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ".bedGraph"
        ),
    shell:
        """
        OUTDIR=$(dirname {output})/alignments
        echo $OUTDIR
        MethylDackel extract {input.genome} {input.alignment} -o $OUTDIR --mergeContext
        """


rule rename_dackel_output:
    input:
        "results/Illumina_pe/{protocol}/result_files/alignments_CpG.bedGraph",
    output:
        "results/Illumina_pe/{protocol}/result_files/methylDackel.bed",
    shell:
        "mv {input} {output}"
