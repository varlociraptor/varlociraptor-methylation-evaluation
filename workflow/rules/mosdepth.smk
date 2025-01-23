rule candidate_splits_to_vcf:
    input:
        "resources/{chro}/candidates_{scatteritem}.bcf",
    output:
        "resources/{chro}/candidates_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/{chro}/candidate_splits_to_vcf_{scatteritem}.log",
    shell:
        """
        bcftools view {input} > {output}
        """


rule candidate_vcf_to_bed:
    input:
        "resources/{chrom}/candidates.vcf",
    output:
        "resources/{chrom}/candidates.bed",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/candidate_vcf_to_bed.py"


# Get coverade over CpG positions
rule mosdepth_bed:
    input:
        bam="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        bai="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        # TODO make chromosome possible as wildcard
        bed=lambda wildcards: expand(
            "resources/{chrom}/candidates.bed",
            chrom=chromosome_by_platform[wildcards.platform],
        ),
    log:
        "logs/mosdepth_bed/{platform}/{protocol}.log",
    output:
        "resources/{platform}/{protocol}/cov.mosdepth.global.dist.txt",
        "resources/{platform}/{protocol}/cov.mosdepth.region.dist.txt",
        "resources/{platform}/{protocol}/cov.regions.bed.gz",
        summary="resources/{platform}/{protocol}/cov.mosdepth.summary.txt",  # this named output is required for prefix parsing
    params:
        extra="--no-per-base --use-median",  # optional
    # additional decompression threads through `--threads`
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v5.5.2/bio/mosdepth"


rule unzip_mosdepth:
    input:
        "resources/{platform}/{protocol}/cov.regions.bed.gz",
    output:
        "resources/{platform}/{protocol}/cov.regions.bed",
    shell:
        """
        gunzip {input}
        """
