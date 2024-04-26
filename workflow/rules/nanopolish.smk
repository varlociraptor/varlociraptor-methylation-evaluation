# Funkitoniert ohne fast5 Files nicht, da wir Nanopolish index aufrufen muessen.


# We need the newer methylation extractor from bsmapz, since the original one is outdated
rule nanopolish_meth_extractor:
    output:
        "resources/ref_tools/nanopolish/calculate_methylation_frequency.py",
    log:
        "../logs/nanopolish_meth_extractor.log",
    params:
        pipeline_path=config["pipeline_path"],
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p mkdir -p $(dirname {output})
        cd $(dirname {output})
        wget https://raw.githubusercontent.com/jts/nanopolish/master/scripts/calculate_methylation_frequency.py
        """


rule compute_nanopore_fastq:
    input:
        "resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
    output:
        "resources/Nanopore/{protocol}/computed_read.fastq",
    shell:
        "samtools fastq {input} > {output}"


rule nanopolish:
    input:
        fastq="resources/Nanopore/{protocol}/computed_read.fastq",
        bam="resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        fasta="resources/genome.fasta",
    output:
        "results/ref_tools/nanopolish/{protocol}/finished.txt",
    params:
        pipeline_path=config["pipeline_path"],
        chromosome=f"chr{str(chromosome_by_platform["Nanopore"])}",
    conda:
        "../envs/nanopolish.yaml"
    shell:
        """
        nanopolish call-methylation -t 8 -r {input.fastq} -b {input.bam} -g {input.fasta} -w "{params.chromosome}" > {params.pipeline_path}results/ref_tools/nanopolish/data1/nanopolish_output.tsv
        touch {output}
        """


rule convert_to_ratio:
    input:
        script="resources/ref_tools/nanopolish/calculate_methylation_frequency.py",
        meth_file="results/ref_tools/nanopolish/{protocol}/finished.txt",
    output:
        "results/ref_tools/nanopolish/{protocol}/methratio.tsv",
    shell:
        """
        {input.script} {input.meth_file} > {output}
        """
