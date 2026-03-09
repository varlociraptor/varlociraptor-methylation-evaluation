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
        alignment="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        jar="resources/ref_tools/Bis-tools/{platform}_{sample}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{platform}_{sample}/genome.fasta",
        genome_index="resources/ref_tools/Bis-tools/{platform}_{sample}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{platform}_{sample}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{platform}_{sample}/alignment.bam.bai",
    log:
        "logs/bissnp/bissnp_prepare/{platform}_{sample}.log",
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


rule bissnp_extract:
    input:
        jar="resources/ref_tools/Bis-tools/{platform}_{sample}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{platform}_{sample}/genome.fasta",
        genome_index="resources/ref_tools/Bis-tools/{platform}_{sample}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{platform}_{sample}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{platform}_{sample}/alignment.bam.bai",
    output:
        cpg=temp("results/single_sample/{platform}/called/{sample}/result_files/cpg_{scatteritem}.raw.vcf"),
        snp=temp("results/single_sample/{platform}/called/{sample}/result_files/snp_{scatteritem}.raw.vcf"),
    conda:
        "../envs/openjdk.yaml"
    params:
        chromosome=chromosome_by_seq_platform.get("Illumina_pe"),
    log:
        "logs/bissnp/bissnp_extract/{platform}_{sample}_{scatteritem}.log",
    benchmark:
        "benchmarks/{platform}/bisSNP/bissnp_extract/{sample}_{scatteritem}.txt"
    threads: 8
    resources:
        mem_mb=64000,
    shell:
        "java -Xmx10G -jar {input.jar} -R {input.genome} -nt {threads} -T BisulfiteGenotyper -I {input.alignment} -vfn1 {output.cpg} -vfn2 {output.snp} -L {params.chromosome} 2> {log}"


rule gather_bisSnp:
    input:
        cpg=gather.split_candidates(
            "results/single_sample/{{platform}}/called/{{sample}}/result_files/cpg_{scatteritem}.raw.vcf",
        ),
        snp=gather.split_candidates(
            "results/single_sample/{{platform}}/called/{{sample}}/result_files/snp_{scatteritem}.raw.vcf",
        ),
    output:
        cpg="results/single_sample/{platform}/called/{sample}/result_files/cpg.raw.vcf",
        snp="results/single_sample/{platform}/called/{sample}/result_files/snp.raw.vcf",
    log:
        "logs/bissnp/gather_bissnp/{platform}_{sample}.log",
    conda:
        "../envs/general.yaml"
    shell:
        """
        cat {input.cpg} > {output.cpg} 2> {log}
        cat {input.snp} > {output.snp} 2>> {log}
        """


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
