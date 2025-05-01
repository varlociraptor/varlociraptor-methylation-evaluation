# Newer versions of BisSNP (0.90.0 or 1.0.0) don't work
# TODO: MissingOutputException after running
rule bissnp_download:
    output:
        "resources/ref_tools/Bis-tools/BisSNP-0.82.2.jar",
    conda:
        "../envs/shell_cmds.yaml"
    log:
        "logs/bissnp/download.log",
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
            chrom=config["platforms"]["Illumina_pe"],
        ),
        genome_index=expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=config["platforms"]["Illumina_pe"],
        ),
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        jar="resources/ref_tools/Bis-tools/{protocol}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{protocol}/genome.fasta",
        genome_index="resources/ref_tools/Bis-tools/{protocol}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{protocol}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{protocol}/alignment.bam.bai",
    log:
        "logs/bissnp/{protocol}/prepare.log",
    shell:
        """
        cp {input.jar} {output.jar} 2> {log}
        cp {input.genome} {output.genome} 2> {log}
        cp {input.genome_index} {output.genome_index} 2> {log}
        cp {input.alignment} {output.alignment} 2> {log}
        cp {input.alignment_index} {output.alignment_index} 2> {log}
        """


rule bissnp_extract_meth:
    input:
        jar="resources/ref_tools/Bis-tools/{protocol}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{protocol}/genome.fasta",
        genome_index="resources/ref_tools/Bis-tools/{protocol}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{protocol}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{protocol}/alignment.bam.bai",
    output:
        cpg="results/Illumina_pe/{protocol}/result_files/cpg.raw.vcf",
        snp="results/Illumina_pe/{protocol}/result_files/snp.raw.vcf",
    conda:
        "../envs/openjdk.yaml"
    params:
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ""
        ),
        chromosome=chromosome_by_platform.get("Illumina_pe"),
    log:
        "logs/bissnp/{protocol}/extract_meth.log",
    shell:
        "java -Xmx10G -jar {input.jar} -R {input.genome} -T BisulfiteGenotyper -I {input.alignment} -vfn1 {output.cpg} -vfn2 {output.snp} -L {params.chromosome} 2> {log}"


# We do not use the official perl script in resources/ref_tools/Bis-tools/utils/vcf2bedGraph.pl because it does not work
# We copied the script and removed line 79: next unless ($splitin[6] eq "PASS" || $splitin[6] eq "Infinity"); because it never triggers
# Furthermore we added coverage information to the output:
#   - line 68: my $head_line = "track type=bedGraph name=${cpg_name_output}.${bissnp_version} description=\"$type methylation level and coverage\" visibility=3";
#   - line 104: my $out_line = "$chr\t$start\t$end\t$methy\t$ct_reads";
rule bissnp_create_bedgraph:
    input:
        perl_script="workflow/scripts/bssnp_bedGraph.pl",
        cpg="results/Illumina_pe/{protocol}/result_files/cpg.raw.vcf",
    output:
        "results/Illumina_pe/{protocol}/result_files/cpg.raw.CG.bedgraph",
    log:
        "logs/bissnp/{protocol}/create_bedgraph.log",
    conda:
        "../envs/openjdk.yaml"
    shell:
        "perl {input.perl_script} {input.cpg} CG 2> {log}"


rule bissnp_rename_output:
    input:
        "results/Illumina_pe/{protocol}/result_files/cpg.raw.CG.bedgraph",
    output:
        "results/Illumina_pe/{protocol}/result_files/bisSNP.bed",
    log:
        "logs/bissnp/{protocol}/rename_output.log",
    shell:
        "mv {input} {output} 2> {log}"
