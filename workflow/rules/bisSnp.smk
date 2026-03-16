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
        output_dir=$(dirname {output})
        mkdir -p "$output_dir"

        # clone only if directory is empty
        if [ -z "$(ls -A "$output_dir")" ]; then
            git clone https://github.com/dnaase/Bis-tools.git "$output_dir"
            wget -O "$output_dir/BisSNP-0.82.2.jar" https://sourceforge.net/projects/bissnp/files/BisSNP-0.82.2/BisSNP-0.82.2.jar/download
        fi

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
        # BisSNP is based on GATK 1.x and requires Bismark-style XR/XG tags to recognize
        # bisulfite reads. bwa-meth BAMs use YC/YD tags instead, which BisSNP does not
        # understand, resulting in 0 callable bases. We therefore use the Bismark-deduplicated
        # BAM as input.
        alignment="resources/ref_tools/bismark/{platform}/dedup/{sample}.deduplicated.bam",
    output:
        jar="resources/ref_tools/Bis-tools/{platform}/{sample}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{platform}/{sample}/genome.fasta",
        genome_index="resources/ref_tools/Bis-tools/{platform}/{sample}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{platform}/{sample}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{platform}/{sample}/alignment.bam.bai",
    log:
        "logs/bissnp/bissnp_prepare/{platform}_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        chromosome=lambda wildcards: chromosome_by_seq_platform.get(wildcards.platform),
    shell:
        """
        cp {input.jar} {output.jar} 2> {log}
        cp {input.genome} {output.genome} 2>> {log}
        cp {input.genome_index} {output.genome_index} 2>> {log}
        samtools sort -o {output.alignment} {input.alignment} 2>> {log}
        samtools index {output.alignment} 2>> {log}
        """


rule bissnp_extract:
    input:
        jar="resources/ref_tools/Bis-tools/{platform}/{sample}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{platform}/{sample}/genome.fasta",
        genome_index="resources/ref_tools/Bis-tools/{platform}/{sample}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{platform}/{sample}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{platform}/{sample}/alignment.bam.bai",
    output:
        cpg=temp("results/single_sample/{platform}/called/{sample}/result_files/cpg.raw.vcf"),
        snp=temp("results/single_sample/{platform}/called/{sample}/result_files/snp.raw.vcf"),
    conda:
        "../envs/openjdk.yaml"
    params:
        chromosome=lambda wildcards: chromosome_by_seq_platform.get(wildcards.platform),
    log:
        "logs/bissnp/bissnp_extract/{platform}_{sample}.log",
    benchmark:
        "benchmarks/{platform}/bisSNP/bissnp_extract/{sample}.txt"
    threads: 8
    resources:
        mem_mb=64000,
    shell:
        "java -Xmx10G -jar {input.jar} -R {input.genome} -nt {threads} -T BisulfiteGenotyper -I {input.alignment} -vfn1 {output.cpg} -vfn2 {output.snp} -L {params.chromosome} 2> {log}"



# We do not use the official perl script in resources/ref_tools/Bis-tools/utils/vcf2bedGraph.pl because it does not work
# We copied the script and removed line 79: next unless ($splitin[6] eq "PASS" || $splitin[6] eq "Infinity"); because it never triggers
# Furthermore we added coverage information to the output:
#   - line 68: my $head_line = "track type=bedGraph name=${cpg_name_output}.${bissnp_version} description=\"$type methylation level and coverage\" visibility=3";
#   - line 104: my $out_line = "$chr\t$start\t$end\t$methy\t$ct_reads";
rule bissnp_create_bedgraph:
    input:
        perl_script=workflow.source_path("../scripts/bissnp_bedGraph.pl"),
        cpg="results/single_sample/{platform}/called/{sample}/result_files/cpg.raw.vcf",
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/cpg.raw.CG.bedgraph",
    log:
        "logs/bissnp/bissnp_create_bedgraph/{platform}_{sample}.log",
    conda:
        "../envs/openjdk.yaml"
    shell:
        "perl {input.perl_script} {input.cpg} CG 2> {log}"


rule bissnp_merge_positions:
    input:
        bedgraph="results/single_sample/{platform}/called/{sample}/result_files/cpg.raw.CG.bedgraph",
        candidates=expand(
            "resources/{chrom}/candidates.bcf",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
        candidates_index=expand(
            "resources/{chrom}/candidates.bcf.csi",
            chrom=config["seq_platforms"].get("Illumina_pe"),
        ),
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/bisSNP.bed",
    log:
        "logs/bissnp/bissnp_merge_positions/{platform}_{sample}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/merge_forward_reverse_positions.py"
