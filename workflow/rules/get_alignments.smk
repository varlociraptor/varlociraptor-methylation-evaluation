chromosome_by_platform = config["platforms"]


rule index_genome:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.bwameth.c2t",
    log:
        "logs/index_genome.log",
    conda:
        "../envs/bwa-meth.yaml"
    shell:
        "bwameth.py index-mem2 {input}"


rule align_reads_pe:
    input:
        fasta_index="resources/genome.fasta.bwameth.c2t",
        fasta="resources/genome.fasta",
        reads1="resources/Illumina_pe/{protocol}/{SRA}/{SRA}_1_trimmed.fastq",
        reads2="resources/Illumina_pe/{protocol}/{SRA}/{SRA}_2_trimmed.fastq",
    output:
        "resources/Illumina_pe/{protocol}/{SRA}/alignment.bam",
    conda:
        "../envs/bwa-meth.yaml"
    log:
        "logs/align_reads_pe_{protocol}_{SRA}.log",
    threads: 30
    resources:
        mem_mb=512,
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    shell:
        """
        bwameth.py --threads {threads} --reference {input.fasta} {input.reads1} {input.reads2}  | samtools view -S -b - > {output}
        """


rule align_reads_se:
    input:
        fasta_index="resources/genome.fasta.bwameth.c2t",
        fasta="resources/genome.fasta",
        reads1="resources/Illumina_se/{protocol}/{SRA}/{SRA}_trimmed.fastq",
    output:
        "resources/Illumina_se/{protocol}/{SRA}/alignment.bam",
    conda:
        "../envs/bwa-meth.yaml"
    log:
        "logs/align_reads_se_{protocol}_{SRA}.log",
    threads: 30
    resources:
        mem_mb=512,
    shell:
        """
        bwameth.py --threads {threads} --reference {input.fasta} | samtools view -S -b - > {output}
        """


rule sort_aligned_reads:
    input:
        "resources/{platform}/{protocol}/{SRA}/alignment.bam",
    output:
        "resources/{platform}/{protocol}/{SRA}/alignment_sorted.bam",
    log:
        "logs/sort_aligned_reads_{platform}_{protocol}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    shell:
        """
        samtools sort -@ {threads}  {input} -o {output}    
        """


rule index_aligned_reads:
    input:
        "resources/{platform}/{protocol}/{SRA}/alignment_sorted.bam",
    output:
        "resources/{platform}/{protocol}/{SRA}/alignment_sorted.bam.bai",
    log:
        "logs/aligned_reads_index_{platform}_{protocol}{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    shell:
        """
        samtools index -@ {threads} {params.pipeline_path}{input}
      """


rule focus_aligned_reads_on_chromosome:
    input:
        bam="resources/{platform}/{protocol}/{SRA}/alignment_sorted.bam",
        index="resources/{platform}/{protocol}/{SRA}/alignment_sorted.bam.bai",
    output:
        bam="resources/{platform}/{protocol}/{SRA}/alignment_focused.bam",
    log:
        "logs/focus_aligned_reads_chrom_{platform}_{protocol}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        chromosome=lambda wildcards: (
            f"chr{chromosome_by_platform[wildcards.platform]}"
            if wildcards.platform == "PacBio" or wildcards.platform == "Nanopore"
            else chromosome_by_platform[wildcards.platform]
        ),
    threads: 10
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    shell:
        """ 
        samtools view -b -o {output.bam} {input} {params.chromosome}
        """


rule filter_mapping_quality:
    input:
        "resources/{platform}/{protocol}/{SRA}/alignment_focused.bam",
    output:
        "resources/{platform}/{protocol}/{SRA}/alignment_focused_filtered.bam",
    log:
        "logs/filter_mapping_quality_{platform}_{protocol}_{SRA}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        min_quality=10,
    threads: 10
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    shell:
        """
        samtools view -q {params.min_quality} -b -o {output} {input}
        """


rule markduplicates_bam:
    input:
        bams="resources/{platform}/{protocol}/{SRA}/alignment_focused_filtered.bam",
    output:
        bam="resources/{platform}/{protocol}/{SRA}/alignment_focused_dedup.bam",
        metrics="resources/{platform}/{protocol}/{SRA}/alignment_focused_dedup.metrics.txt",
    log:
        "logs/markduplicates_bam__{platform}_{protocol}_{SRA}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024,
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    wrapper:
        "v2.6.0/bio/picard/markduplicates"


rule merge_bams:
    input:
        get_protocol_sra,
    output:
        "resources/{platform, [^/]+}/{protocol,[^/]+}/alignment_focused_dedup.bam",
    log:
        "logs/merge_bams_{platform}_{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    shell:
        """
        echo {input}
        samtools merge {output} {input}
        """


rule downsample_bams:
    input:
        "resources/{platform}/{protocol}/alignment_focused_dedup.bam",
    output:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup.bam",
    log:
        "logs/downsample_bams_{platform}_{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    shell:
        """
        samtools view -s 0.99 -b -o {output} {input}
        """


rule aligned_downsampled_index:
    input:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup.bam",
    output:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup.bam.bai",
    log:
        "logs/aligned_downsampled_index_{platform}_{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    shell:
        """
        samtools index -@ {threads} {params.pipeline_path}{input}
        """


# Die Regel funktioniert nur ausserhalb von Snakemake?
rule bam_to_sam:
    input:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup.bam",
    output:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup.sam",
    log:
        "logs/bam_to_sam_{platform}_{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    shell:
        """
        samtools view -@ {threads} -h  -o {params.pipeline_path}/{output}  {params.pipeline_path}{input}   

        """


rule rename_chromosomes_in_sam:
    input:
        sam="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup.sam",
    output:
        renamed="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.sam",
    log:
        "logs/rename_chromosomes_in_sam_{platform}_{protocol}.log",
    conda:
        "../envs/python.yaml"
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    script:
        "../scripts/rename_alignment.py"


rule sam_to_bam:
    input:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.sam",
    output:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
    log:
        "logs/sam_to_bam_{platform}_{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    wildcard_constraints:
        protocol="^(?!simulated_data$).*",
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """
        samtools view -bS {input} > {output}   
        """


rule aligned_downsampled_reads_dedup_index:
    input:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
    output:
        "resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    log:
        "logs/aligned_downsampled_reads_dedup_index_{platform}_{protocol}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        samtools index -@ {threads} {params.pipeline_path}{input}
        """


rule aligned_reads_candidates_region:
    input:
        alignment="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        index="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        candidate=lambda wildcards: expand(
            "resources/{chrom}/candidates_{{scatteritem}}.bcf",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
    output:
        "resources/{platform}/{protocol}/candidate_specific/alignment_{scatteritem}.bam",
    log:
        "logs/aligned_reads_candidates_region_{platform}_{protocol}_{scatteritem}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        window_size=config["max_read_length"],
        chromosome=lambda wildcards: chromosome_by_platform[wildcards.platform],
    shell:
        """
        set +o pipefail;
        start=$(bcftools query -f '%POS\\n' {input.candidate} | head -n 1)
        end=$(bcftools query -f '%POS\\n' {input.candidate} | tail -n 1)
        end=$((end + {params.window_size}))
        samtools view -b {input.alignment} "{params.chromosome}:$start-$end" > {output}
        """


# TODO: This rule is actually unnecessary and consumes a lot of time. The problem is that Varlociraptor does not work on bam files that have no reads. Due to the previous step, it can happen that all candidates are outside the reads and the bam file therefore has no reads. Therefore we simply add the last read of the original bam file for each bam file.
rule aligned_reads_candidates_region_valid:
    input:
        alignment_orig="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_specific="resources/{platform}/{protocol}/candidate_specific/alignment_{scatteritem}.bam",
    output:
        bam="resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam",
    log:
        "logs/aligned_reads_candidates_region_valid_{platform}_{protocol}_{scatteritem}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        last_read=$(samtools view {input.alignment_orig} | tail -n 1)
        echo -e "$(samtools view -H {input.alignment_specific})\n$last_read" | samtools view -bS - > {output.bam}
        """


rule aligned_reads_candidates_region_index:
    input:
        "resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam",
    output:
        "resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam.bai",
    log:
        "logs/aligned_reads_candidates_region_index_{platform}_{protocol}_{scatteritem}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        samtools index -@ {threads} {params.pipeline_path}{input}
        """
