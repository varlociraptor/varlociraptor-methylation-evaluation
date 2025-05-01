rule mason_download:
    output:
        mason_dir=directory("resources/tools/seqan/apps/mason2"),
        mason="resources/tools/seqan/apps/mason2/methylation_levels.h",
    log:
        "../logs/mason/download.log",
    conda:
        "../envs/shell_cmds.yaml"
    shell:
        """
        mkdir -p resources/tools 2> {log}
        cd resources/tools 2> {log}
        git clone git@github.com:seqan/seqan.git 2> {log}
        """


rule mason_fake_methylation:
    input:
        chrom="resources/chromosome_{chrom}.fasta",
    output:
        methylation="resources/Illumina_pe/simulated_data/chromosome_{chrom}_meth.fa",
    conda:
        "../envs/mason.yaml"
    log:
        "logs/mason_methylation/fake_methylation_{chrom}.log",
    shell:
        """
        mason_methylation --in {input.chrom} \
            --methylation-levels \
            --meth-cg-sigma 0.3 \
            --meth-cg-mu 0.5 \
            --out {output.methylation}  2> {log}
        """


rule mason_fake_variants:
    input:
        "resources/chromosome_{chrom}.fasta",
    output:
        "resources/Illumina_pe/simulated_data/chromosome_{chrom}_variants.vcf",
    conda:
        "../envs/mason.yaml"
    log:
        "logs/mason_variants/fake_variants_{chrom}.log",
    shell:
        """
        mason_variator --in-reference {input} \
            --out-vcf {output}  2> {log}
        """
        # --snp-rate 0.01 \


rule mason_fake_reads:
    input:
        genome="resources/chromosome_{chrom}.fasta",
        variants="resources/Illumina_pe/simulated_data/chromosome_{chrom}_variants.vcf",
        methylation="resources/Illumina_pe/simulated_data/chromosome_{chrom}_meth.fa",
    output:
        f1="resources/Illumina_pe/simulated_data/chromosome_{chrom}_f1.fastq",
        f2="resources/Illumina_pe/simulated_data/chromosome_{chrom}_f2.fastq",
    conda:
        "../envs/mason.yaml"
    log:
        "logs/mason_reads/fake_reads_{chrom}.log",
    params:
        num_fragments=config["num_simulated_reads"],
    shell:
        """
        mason_simulator --input-reference {input.genome} \
                --input-vcf {input.variants} \
                --num-fragments {params.num_fragments} \
                --out {output.f1} \
                --out-right {output.f2} \
                --meth-fasta-in {input.methylation} \
                --enable-bs-seq \
                --illumina-read-length 150  2> {log}
        """


rule mason_align_reads:
    input:
        fasta=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        fasta_index=expand(
            "resources/chromosome_{chrom}.fasta.bwameth.c2t",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        f1=expand(
            "resources/Illumina_pe/simulated_data/chromosome_{chrom}_f1.fastq",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        f2=expand(
            "resources/Illumina_pe/simulated_data/chromosome_{chrom}_f2.fastq",
            chrom=config["platforms"]["Illumina_pe"],
        ),
    output:
        "resources/Illumina_pe/simulated_data/alignment.sam",
    conda:
        "../envs/bwa-meth.yaml"
    log:
        "logs/mason/align_reads.log",
    threads: 30
    shell:
        """
        set -o pipefail
        bwameth.py index-mem2 {input.fasta} 2> {log}
        bwameth.py --threads {threads} --reference {input.fasta} {input.f1} {input.f2} > {output} 2> {log}
        """


rule mason_sam_to_bam:
    input:
        "resources/Illumina_pe/simulated_data/alignment.sam",
    output:
        "resources/Illumina_pe/simulated_data/alignment.bam",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/mason/sam_to_bam.log",
    threads: 30
    shell:
        "samtools view -Sb {input} > {output} 2> {log}"


rule mason_sort_reads:
    input:
        "resources/Illumina_pe/simulated_data/alignment.bam",
    output:
        # Name it like that in order to skip filtering on qual, mark_duplicates, ...
        "resources/Illumina_pe/simulated_data/alignment_focused_downsampled_dedup_renamed.bam",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/mason/sort_reads.log",
    threads: 10
    shell:
        "samtools sort -@ {threads}  {input} -o {output} 2> {log}"


# Mason has a different meth ratio for forward and reverse strands.
# That is why we need to compute the coverage on the forward and reverse strand independently.
rule mason_alignment_forward:
    input:
        "resources/Illumina_pe/simulated_data/alignment.bam",
    output:
        first="resources/Illumina_pe/simulated_data/alignment_99.bam",
        second="resources/Illumina_pe/simulated_data/alignment_147.bam",
        forward="resources/Illumina_pe/simulated_data/alignment_forward.bam",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/mason/alignment_forward.log",
    shell:
        """
        samtools view -b -f 64 -F 16 {input} > {output.first} 2> {log}
        samtools view -b -f 16 -F 64 {input} > {output.second} 2> {log}
        samtools merge {output.forward} {output.first} {output.second} 2> {log}
        """


rule mason_alignment_reverse:
    input:
        "resources/Illumina_pe/simulated_data/alignment.bam",
    output:
        first="resources/Illumina_pe/simulated_data/alignment_83.bam",
        second="resources/Illumina_pe/simulated_data/alignment_163.bam",
        rev="resources/Illumina_pe/simulated_data/alignment_reverse.bam",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/mason/alignment_reverse.log",
    shell:
        """
        samtools view -b -f 16 -f 64 {input} > {output.first} 2> {log}
        samtools view -b -F 16 -F 64 {input} > {output.second} 2> {log}
        samtools merge {output.rev} {output.first} {output.second} 2> {log}
        """


rule mason_sort_oriented_reads:
    input:
        "resources/Illumina_pe/simulated_data/alignment_{orientation}.bam",
    output:
        "resources/Illumina_pe/simulated_data/alignment_sorted_{orientation}.bam",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    log:
        "logs/mason/sort_oriented_reads_{orientation}.log",
    shell:
        "samtools sort -@ {threads}  {input} -o {output} 2> {log}"


rule mason_index_oriented_alignment:
    input:
        "resources/Illumina_pe/simulated_data/alignment_sorted_{orientation}.bam",
    output:
        "resources/Illumina_pe/simulated_data/alignment_sorted_{orientation}.bam.bai",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/mason/index_oriented_alignment_{orientation}.log",
    shell:
        "samtools index {input} 2> {log}"


rule mason_coverage:
    input:
        bam="resources/Illumina_pe/simulated_data/alignment_sorted_{orientation}.bam",
        bai="resources/Illumina_pe/simulated_data/alignment_sorted_{orientation}.bam.bai",
        bed=expand(
            "resources/{chrom}/candidates.bed",
            chrom=config["platforms"]["Illumina_pe"],
        ),
    output:
        "resources/Illumina_pe/simulated_data/{orientation}_cov.mosdepth.global.dist.txt",
        "resources/Illumina_pe/simulated_data/{orientation}_cov.mosdepth.region.dist.txt",
        "resources/Illumina_pe/simulated_data/{orientation}_cov.regions.bed.gz",
        summary="resources/Illumina_pe/simulated_data/{orientation}_cov.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mason/coverage_{orientation}.log",
    params:
        extra="--no-per-base --use-median",  # optional
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v5.5.2/bio/mosdepth"


rule mason_unzip_coverage:
    input:
        "resources/Illumina_pe/simulated_data/{orientation}_cov.regions.bed.gz",
    output:
        "resources/Illumina_pe/simulated_data/{orientation}_cov.regions.bed",
    log:
        "logs/mason/unzip_coverage_{orientation}.log",
    shell:
        "gunzip {input} 2> {log}"


rule mason_compute_truth:
    input:
        cov_forward="resources/Illumina_pe/simulated_data/forward_cov.regions.bed",
        cov_reverse="resources/Illumina_pe/simulated_data/reverse_cov.regions.bed",
        methylation="resources/Illumina_pe/simulated_data/chromosome_{chrom}_meth.fa",
        candidates="resources/{chrom}/candidates.vcf",
    output:
        "resources/Illumina_pe/simulated_data/chromosome_{chrom}_truth.bed",
    conda:
        "../envs/python.yaml"
    log:
        "logs/mason/compute_truth_{chrom}.log",
    script:
        "../scripts/mason_ascii_to_meth.py"
