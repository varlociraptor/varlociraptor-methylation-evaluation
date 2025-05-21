# Computes the coverage over the input alignment in order to have a common coverage over which we can stratidfy the plots of the different tools
rule candidate_to_vcf:
    input:
        "resources/{chrom}/candidates_{scatteritem}.bcf",
    output:
        "resources/{chrom}/candidates_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/mosdepth/{chrom}/candidate_to_vcf_{scatteritem}.log",
    shell:
        "bcftools view {input} > {output} 2> {log}"


rule vcf_to_bed:
    input:
        "resources/{chrom}/candidates.vcf",
    output:
        "resources/{chrom}/candidates.bed",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/mosdepth/{chrom}/vcf_to_bed.log",
    script:
        "../scripts/vcf_to_bed.py"


# Get coverade over CpG positions
rule mosdepth_compute_coverage:
    input:
        bam="resources/{seq_platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        bai="resources/{seq_platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        bed=lambda wildcards: expand(
            "resources/{chrom}/candidates.bed",
            chrom=chromosome_by_seq_platform[wildcards.seq_platform],
        ),
    log:
        "logs/mosdepth/{seq_platform}/{protocol}/compute_coverage.log",
    output:
        "resources/{seq_platform}/{protocol}/cov.mosdepth.global.dist.txt",
        "resources/{seq_platform}/{protocol}/cov.mosdepth.region.dist.txt",
        "resources/{seq_platform}/{protocol}/cov.regions.bed.gz",
        summary="resources/{seq_platform}/{protocol}/cov.mosdepth.summary.txt",  # this named output is required for prefix parsing
    params:
        extra="--no-per-base --use-median",  # optional
    # additional decompression threads through `--threads`
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v5.5.2/bio/mosdepth"


# We can't just do gunzip because of symlinks from the cluster
rule mosdepth_unzip_results:
    input:
        "resources/{seq_platform}/{protocol}/cov.regions.bed.gz",
    output:
        "resources/{seq_platform}/{protocol}/cov.regions.bed",
    log:
        "logs/mosdepth/{seq_platform}/{protocol}/unzip_results.log",
    shell:
        """
        cp {input} tmpfile.gz
        gunzip -k tmpfile.gz
        mv tmpfile {output}
        """
