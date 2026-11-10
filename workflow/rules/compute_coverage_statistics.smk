# Our candidate file spans CG positions, we are only interested on the coverage of C positions
rule focus_candidates_on_c:
    input:
        "resources/{chrom}/candidates.bed",
    output:
        "resources/{chrom}/candidates_focus_c.bed",
    log:
        "logs/mason/focus_candidates_on_c/{chrom}.log",
    conda:
        "../envs/general.yaml"
    shell:
        r"""
        awk 'BEGIN{{OFS="\t"}}{{
            center=$2+1
            print $1, center, center+1
        }}' {input} > {output}
        """


rule compute_coverage:
    input:
        bam="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        bai="resources/{seq_platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        bed=lambda wildcards: expand(
            "resources/{chrom}/candidates_focus_c.bed",
            chrom=config["seq_platforms"].get(wildcards.seq_platform, []),
        ),
    output:
        "results/{call_type}/{seq_platform}/coverages/{sample}_{mapq}.mosdepth.global.dist.txt",
        "results/{call_type}/{seq_platform}/coverages/{sample}_{mapq}.mosdepth.region.dist.txt",
        "results/{call_type}/{seq_platform}/coverages/{sample}_{mapq}.regions.bed.gz",
        summary="results/{call_type}/{seq_platform}/coverages/{sample}_{mapq}.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mason/mason_coverage/{call_type}_{seq_platform}_{sample}_{mapq}.log",
    wildcard_constraints:
        sample="(?!all_samples).*",
    threads: 4  # This value - 1 will be sent to `--threads`
    params:
        extra=lambda wildcards: (
            f"--no-per-base --use-median --mapq {wildcards.mapq}"
            if wildcards.mapq != "all"
            else ""
        ),
        # optional
    wrapper:
        "v5.5.2/bio/mosdepth"


rule common_coverage_Illumina:
    input:
        lambda wildcards: expand(
            "results/single_sample/Illumina_pe/coverages/{sample}_{{REP}}_{{mapq}}.regions.bed.gz",
            sample=config["samples"]["Illumina_pe"],
        ),
    output:
        "results/single_sample/Illumina_pe/coverages/all_samples_{REP}_{mapq}.regions.bed.gz",
    log:
        "logs/common_coverage_Illumina/all_samples_{REP}_{mapq}.log",
    conda:
        "../envs/python.yaml"
    threads: 1
    script:
        "../scripts/merge_coverage.py"
