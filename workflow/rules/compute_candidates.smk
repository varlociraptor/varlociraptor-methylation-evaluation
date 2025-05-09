rule find_candidates:
    input:
        varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
        fasta="resources/chromosome_{chromosome}.fasta",
    output:
        "resources/{chromosome}/candidates.bcf",
    log:
        "logs/candidates/find_candidates_{chromosome}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} methylation-candidates {input.fasta} {output} 2> {log}"


# Use only those candidates, which are in the bam file
rule relevant_positions:
    input:
        bam="resources/{seq_platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        index="resources/{seq_platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        bcf="resources/{chrom}/candidates.bcf",
    output:
        "resources/{seq_platform}/{protocol}/{chrom}/candidates_shortened.bcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/candidates/{seq_platform}/{protocol}/relevant_positions_{chrom}.log",
    shell:
        """
        set +o pipefail;
        start=$(samtools depth {input.bam} | awk '$3 > 0 {{print $2}}' | head -n 1) 2> {log}
        end=$(samtools depth {input.bam} | awk '$3 > 0 {{print $2}}' | tail -n 1) 2> {log}
        mkdir -p $(dirname {output}) 2> {log}
        bcftools index -f {input.bcf} 2> {log}
        bcftools view -r $(samtools idxstats {input.bam} | awk '$3 > 0 {{print $1}}'):${{start}}-${{end}} {input.bcf} -o {output} 2> {log}
        """


rule split_candidates:
    input:
        "resources/{seq_platform}/{protocol}/{chrom}/candidates_shortened.bcf",
    output:
        scatter.split_candidates(
            "resources/{{seq_platform}}/{{protocol}}/{{chrom}}/candidates_{scatteritem}.bcf"
        ),
    log:
        "logs/candidates/{seq_platform}/{protocol}/split_candidates_{chrom}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output} 2> {log}"


# This is only to debug the adjusted mapq value. Will stay until completely clarified
# rule debug_filter_candidates_on_mapq:
#     input:
#         alignment="resources/{seq_platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam",
#         index="resources/{seq_platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam.bai",
#         candidates="resources/{seq_platform}/{protocol}/{chrom}/candidates_{scatteritem}.bcf",
#     output:
#         "resources/{seq_platform}/{protocol}/{chrom}/candidates_filtered_mapq_{scatteritem}.bcf",
#     log:
#         "logs/filter_candidates_on_mapq_{seq_platform}_{protocol}_{chrom}_{scatteritem}.log",
#     conda:
#         "../envs/plot.yaml"
#     script:
#         "../scripts/debug_filter_candidates_on_mapq.py"
