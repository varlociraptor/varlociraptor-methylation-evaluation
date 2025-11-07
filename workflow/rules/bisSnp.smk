# Newer versions of BisSNP (0.90.0 or 1.0.0) don't work
# TODO: MissingOutputException after running on cluster
rule bissnp_download:
    output:
        "resources/ref_tools/Bis-tools/BisSNP-0.82.2.jar",
    conda:
        "../envs/shell_cmds.yaml"
    log:
        "logs/bissnp/bissnp_download/download.log",
    shell:
        """
        touch {log}
        mkdir -p resources/ref_tools
        cd resources/ref_tools
        git clone https://github.com/dnaase/Bis-tools.git
        cd Bis-tools
        wget -O BisSNP-0.82.2.jar https://sourceforge.net/projects/bissnp/files/BisSNP-0.82.2/BisSNP-0.82.2.jar/download
        """


# All files need to be in the same dir
rule bissnp_prepare:
    input:
        jar="resources/ref_tools/Bis-tools/BisSNP-0.82.2.jar",
        genome=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        genome_index=expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        alignment="resources/Illumina_pe/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        jar="resources/ref_tools/Bis-tools/{sample}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{sample}/genome.fasta",
        genome_index="resources/ref_tools/Bis-tools/{sample}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{sample}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{sample}/alignment.bam.bai",
    log:
        "logs/bissnp/bissnp_prepare/{sample}.log",
    conda:
        "../envs/general.yaml"
    shell:
        """
        cp {input.jar} {output.jar} 2> {log}
        cp {input.genome} {output.genome} 2> {log}
        cp {input.genome_index} {output.genome_index} 2> {log}
        cp {input.alignment} {output.alignment} 2> {log}
        cp {input.alignment_index} {output.alignment_index} 2> {log}
        """


# We tried to split the alignment file in smaller parts to process faster but it did not work properly
# rule split_bisSNP_alignments:
#     input:
#         alignment="resources/ref_tools/Bis-tools/{sample}/alignment.bam",
#         alignment_index="resources/ref_tools/Bis-tools/{sample}/alignment.bam.bai",
#     output:
#         bams=expand(
#             "resources/ref_tools/Bis-tools/{{sample}}/alignment_{i}-of-{n}.bam",
#             i=range(1, config["scatter_number"] + 1),
#             n=config["scatter_number"],
#         ),
#         header="resources/ref_tools/Bis-tools/{sample}/header.sam",
#     log:
#         "logs/bissnp/split_bisSNP_alignments/{sample}.log",
#     conda:
#         "../envs/samtools.yaml"
#     params:
#         n=lambda wildcards: config["scatter_number"],
#     shell:
#         r"""
#         total=$(samtools view -c {input.alignment})
#         per_part=$(( (total + {params.n} - 1) / {params.n} ))

#         echo "Total reads: $total" > {log}
#         echo "Splitting into {params.n} parts (~$per_part reads each)" >> {log}

#         if [ {params.n} -eq 1 ]; then
#             echo "Only one part requested — copying input BAM directly." >> {log}
#             cp {input.alignment} {output[0]}
#             samtools index {output[0]}
#             echo "Done." >> {log}
#             exit 0
#         fi

#         samtools view -H {input.alignment} > header.sam
#         samtools view {input.alignment} | \
#             awk -v per_part=$per_part -v prefix="resources/ref_tools/Bis-tools/{wildcards.sample}/alignment_" -v n={params.n} \
#             'BEGIN {{file_index=1; line_count=0}}
#              {{line_count++; print >> (prefix file_index "-of-" n ".sam");
#               if (line_count >= per_part) {{close(prefix file_index "-of-" n ".sam"); file_index++; line_count=0}}}}
#              END {{for (i=file_index; i<=n; i++) close(prefix i "-of-" n ".sam")}}'

#         for i in $(seq 1 {params.n}); do
#             cat header.sam resources/ref_tools/Bis-tools/{wildcards.sample}/alignment_${{i}}-of-{params.n}.sam | \
#                 samtools view -b -o resources/ref_tools/Bis-tools/{wildcards.sample}/alignment_${{i}}-of-{params.n}.bam -
#             rm resources/ref_tools/Bis-tools/{wildcards.sample}/alignment_${{i}}-of-{params.n}.sam
#             samtools index resources/ref_tools/Bis-tools/{wildcards.sample}/alignment_${{i}}-of-{params.n}.bam
#         done

#         rm header.sam
#         echo "Splitting into {params.n} BAMs completed." >> {log}
#         """


# rule index_bisSNP_alignments:
#     input:
#         "resources/ref_tools/Bis-tools/{sample}/alignment.bam",
#     output:
#         "resources/ref_tools/Bis-tools/{sample}/alignment.bam.bai",
#     log:
#         "logs/bissnp/index_bisSNP_alignments/{sample}.log",
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         """
#         samtools index {input} 2> {log}
#         """


rule bissnp_extract:
    input:
        jar="resources/ref_tools/Bis-tools/{sample}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{sample}/genome.fasta",
        genome_index="resources/ref_tools/Bis-tools/{sample}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{sample}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{sample}/alignment.bam.bai",
    output:
        cpg="results/single_sample/Illumina_pe/called/{sample}/result_files/cpg.raw.vcf",
        snp="results/single_sample/Illumina_pe/called/{sample}/result_files/snp.raw.vcf",
    conda:
        "../envs/openjdk.yaml"
    params:
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ""
        ),
        chromosome=chromosome_by_seq_platform.get("Illumina_pe"),
    log:
        "logs/bissnp/bissnp_extract/{sample}.log",
    benchmark:
        "benchmarks/Illumina_pe/bisSNP/bissnp_extract/{sample}.txt"
    resources:
        mem_mb=64000,
    shell:
        "java -Xmx10G -jar {input.jar} -R {input.genome} -T BisulfiteGenotyper -I {input.alignment} -vfn1 {output.cpg} -vfn2 {output.snp} -L {params.chromosome} 2> {log}"


rule gather_bisSnp:
    input:
        cpg=gather.split_candidates(
            "results/single_sample/Illumina_pe/called/{{sample}}/result_files/cpg.raw.vcf",
        ),
        snp=gather.split_candidates(
            "results/single_sample/Illumina_pe/called/{{sample}}/result_files/snp.raw.vcf",
        ),
    output:
        cpg="results/single_sample/Illumina_pe/called/{sample}/result_files/cpg.raw.vcf",
        snp="results/single_sample/Illumina_pe/called/{sample}/result_files/snp.raw.vcf",
    log:
        "logs/bissnp/gather_bissnp/{sample}.log",
    conda:
        "../envs/general.yaml"
    shell:
        """
        cat {input.cpg} > {output.cpg} 2> {log}
        cat {input.snp} > {output.snp} 2> {log}
        """


# We do not use the official perl script in resources/ref_tools/Bis-tools/utils/vcf2bedGraph.pl because it does not work
# We copied the script and removed line 79: next unless ($splitin[6] eq "PASS" || $splitin[6] eq "Infinity"); because it never triggers
# Furthermore we added coverage information to the output:
#   - line 68: my $head_line = "track type=bedGraph name=${cpg_name_output}.${bissnp_version} description=\"$type methylation level and coverage\" visibility=3";
#   - line 104: my $out_line = "$chr\t$start\t$end\t$methy\t$ct_reads";
rule bissnp_create_bedgraph:
    input:
        perl_script="../scripts/bissnp_bedGraph.pl",
        cpg="results/single_sample/Illumina_pe/called/{sample}/result_files/cpg.raw.vcf",
    output:
        "results/single_sample/Illumina_pe/called/{sample}/result_files/cpg.raw.CG.bedgraph",
    log:
        "logs/bissnp/bissnp_create_bedgraph/{sample}.log",
    conda:
        "../envs/openjdk.yaml"
    shell:
        "perl {input.perl_script} {input.cpg} CG 2> {log}"


rule bissnp_merge_positions:
    input:
        bedgraph="results/single_sample/Illumina_pe/called/{sample}/result_files/cpg.raw.CG.bedgraph",
        candidates=expand(
            "resources/{chrom}/candidates.bcf",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        candidates_index=expand(
            "resources/{chrom}/candidates.bcf.csi",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
    output:
        "results/single_sample/Illumina_pe/called/{sample}/result_files/bisSNP.bed",
    log:
        "logs/bissnp/bissnp_merge_positions/{sample}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/merge_forward_reverse_positions.py"
