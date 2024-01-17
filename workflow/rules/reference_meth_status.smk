rule methylDackel:
    input:
        chromosome=expand("resources/chromosome_{chromosome}.fasta", chromosome=chromosome_conf["chromosome"]),
        alignment="resources/Illumina_{read_type}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_{read_type}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        "results/Illumina_{read_type}/{protocol}/alignments_CpG.bedGraph",
    conda:
        "../envs/methylDackel.yaml"
    log:
        "logs/methylDackel_{read_type}_{protocol}.log",
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """
        MethylDackel extract {input.chromosome} {input.alignment} -o {params.pipeline_path}/results/Illumina_{wildcards.read_type}/{wildcards.protocol}/alignments --mergeContext
        """

rule pb_CpG_tools:
    input:
        pb_tools_dir=directory("../pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu"),
        alignment="resources/PacBio/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/PacBio/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        "results/PacBio/{protocol}/alignments_CpG.combined.bed"
    log:
        "logs/pb_CpG_tools_{protocol}.log",
    params:
        base_dir=config["base_dir"],
        prefix="results/PacBio/{protocol}/alignments_CpG"
    shell:
        """
        {params.base_dir}pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
        --bam {input.alignment} \
        --output-prefix {params.prefix} \
        --model {params.base_dir}pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
        --threads 8
        """

# Needs a fasta with >chr1 instead of >1
rule modkit:
    input:
        installed="modkit_installed.flag",
        alignment="resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup.bam",
        chromosome=expand("resources/chromosome_{chromosome}.fasta", chromosome=chromosome_conf["chromosome"])
    output:
        "results/Nanopore/{protocol}/alignments_CpG.combined.bed"
    conda:
        "../envs/modkit.yaml"
    log:
        "logs/modkit_{protocol}.log",
    shell:
        """
        export PATH=$PATH:~/.cargo/bin
        export PATH=$PATH:/homes/aprinz/.cargo/bin
        modkit pileup {input.alignment} {output} --cpg --ref {input.chromosome} --force-allow-implicit --combine-strands --log-filepath log_modkit.txt
        """




# rule nanopore_bam_to_fastq:
#     input:
#         "resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
#     output:
#         fastq1="resources/Nanopore/{protocol}_R1.fastq",
#         fastq2="resources/Nanopore/{protocol}_R2.fastq",
#     log:
#         "logs/picard/sam_to_fastq/{protocol}.log",
#     params:
#         extra="",  # optional: Extra arguments for picard.
#     # optional specification of memory usage of the JVM that snakemake will respect with global
#     # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
#     # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
#     # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
#     resources:
#         mem_mb=1024,
#     wrapper:
#         "v3.3.3/bio/picard/samtofastq"


# rule combine_nanopore_fastqs:
#     input:        
#         fastq1="resources/Nanopore/{protocol}_R1.fastq",
#         fastq2="resources/Nanopore/{protocol}_R2.fastq",
#     output:
#         "resources/Nanopore/{protocol}_combined.fastq",
#     log:
#         "logs/combine_nanopore_fastqs_{protocol}"
#     shell:
#         "cat {input.fastq1} {input.fastq2} > {output}"
    
# rule call_nanopolish:
#     input:
#         fastq="resources/Nanopore/{protocol}_combined.fastq",
#         alignment="resources/Nanopore/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
#         genome="resources/chromosome.fasta"
#     output:
#         "results/Nanopore/{protocol}/alignments_CpG.tsv"
#     conda:
#         "envs/nanopolish.yaml"
#     log:
#         "logs/call_nanopolish_{protocol}.log"
#     params:
#         threads=8,
#         chr=chromosome_conf["chromosome"],
#     shell:
#         """
#         nanopolish call-methylation -t 8 -r {input.fastq} -b {input.alignment} -g {input.genome} -w "{params.chr}" > {output}
#         """

#         # nanopolish index -d [Pfad_zu_den_Nanopore-Rohdaten] resources/Nanopore/data1_combined.fastq

