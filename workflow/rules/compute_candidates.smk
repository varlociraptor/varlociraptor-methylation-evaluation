rule find_candidates:
    input:
        varlo=directory("resources/tools/varlociraptor"),
        fasta="resources/chromosome_{chromosome}.fasta",
    output:
        "resources/{chromosome}/candidates.bcf",
    log:
        "logs/{chromosome}/find_candidates.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        """ 
        cd {input.varlo}
        cargo run -- methylation-candidates {input.fasta} {output}
        """


# Use only those candidates, which are in the bam file
rule relevant_positions:
    input:
        bam="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        index="resources/{platform}/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        bcf="resources/{chrom}/candidates.bcf",
    output:
        "resources/{platform}/{protocol}/{chrom}/candidates_shortened.bcf",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        set +o pipefail;
        start=$(samtools depth {input.bam} | awk '$3 > 0 {{print $2}}' | head -n 1)
        end=$(samtools depth {input.bam} | awk '$3 > 0 {{print $2}}' | tail -n 1)
        mkdir -p $(dirname {output})
        bcftools index -f {input.bcf}
        bcftools view -r $(samtools idxstats {input.bam} | awk '$3 > 0 {{print $1}}'):${{start}}-${{end}} {input.bcf} -o {output}
        """


rule split_candidates:
    input:
        "resources/{platform}/{protocol}/{chrom}/candidates_shortened.bcf",
    output:
        scatter.split_candidates(
            "resources/{{platform}}/{{protocol}}/{{chrom}}/candidates_{scatteritem}.bcf"
        ),
    # log:
    #     "logs/{chrom}split_candidates.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"


# This is only to debug the adjusted mapq value. Will stay until completely clarified
# rule filter_candidates_on_mapq:
#     input:
#         alignment="resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam",
#         index="resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam.bai",
#         candidates="resources/{platform}/{protocol}/{chrom}/candidates_{scatteritem}.bcf",
#     output:
#         "resources/{platform}/{protocol}/{chrom}/candidates_filtered_mapq_{scatteritem}.bcf",
#     log:
#         "logs/filter_candidates_on_mapq_{platform}_{protocol}_{chrom}_{scatteritem}.log",
#     conda:
#         "../envs/plot.yaml"
#     script:
#         "../scripts/filter_candidates_on_mapq.py"
