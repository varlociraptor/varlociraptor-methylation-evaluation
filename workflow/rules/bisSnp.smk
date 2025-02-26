# Newer versions of BisSNP (0.90.0 or 1.0.0) don't work
rule download_BisSNP:
    output:
        "resources/ref_tools/Bis-tools/BisSNP-0.82.2.jar",
    log:
        "../logs/download_BisSNP.log",
    params:
        pipeline_path=config["pipeline_path"],
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p resources/ref_tools
        cd /resources/ref_tools
        git clone https://github.com/dnaase/Bis-tools.git
        cd Bis-tools
        wget -O BisSNP-0.82.2.jar https://sourceforge.net/projects/bissnp/files/BisSNP-0.82.2/BisSNP-0.82.2.jar/download
        """


# All files need to be in the same dir
rule prepare_bissnp:
    input:
        jar="resources/ref_tools/Bis-tools/BisSNP-0.82.2.jar",
        genome="resources/genome.fasta",
        genome_index="resources/genome.fasta.fai",
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup.bam.bai",
    output:
        jar="resources/ref_tools/Bis-tools/{protocol}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{protocol}/genome.fasta",
        genome_index="resources/ref_tools/Bis-tools/{protocol}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{protocol}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{protocol}/alignment.bam.bai",
    log:
        "logs/prepare_bissnp_{protocol}.log",
    shell:
        """
        cp {input.jar} {output.jar}
        cp {input.genome} {output.genome}
        cp {input.genome_index} {output.genome_index}
        cp {input.alignment} {output.alignment}
        cp {input.alignment_index} {output.alignment_index}
        """


rule BisSNP_extract:
    input:
        jar="resources/ref_tools/Bis-tools/{protocol}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{protocol}/genome.fasta",
        genome_index="resources/ref_tools/Bis-tools/{protocol}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{protocol}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{protocol}/alignment.bam.bai",
    output:
        cpg="results/Illumina_pe/{protocol}/result_files/cpg.raw.vcf",
        snp="results/Illumina_pe/{protocol}/result_files/snp.raw.vcf",
    log:
        "logs/pb_CpG_tools_{protocol}.log",
    conda:
        "../envs/bisSnp.yaml"
    params:
        base_dir=config["base_dir"],
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ""
        ),
        chromosome=chromosome_by_platform.get("Illumina_pe"),
    shell:
        """
        java -Xmx10G -jar {input.jar} -R {input.genome} -T BisulfiteGenotyper -I {input.alignment} -vfn1 {output.cpg} -vfn2 {output.snp} -L {params.chromosome}
        """


# We do not use the official perl script in resources/ref_tools/Bis-tools/utils/vcf2bedGraph.pl because it does not work
# We copied the script and removed line 79: next unless ($splitin[6] eq "PASS" || $splitin[6] eq "Infinity"); because it never triggers
# Furthermore we added coverage information to the output:
#   - line 68: my $head_line = "track type=bedGraph name=${cpg_name_output}.${bissnp_version} description=\"$type methylation level and coverage\" visibility=3";
#   - line 104: my $out_line = "$chr\t$start\t$end\t$methy\t$ct_reads";
rule BisSNP_bedGraph:
    input:
        perl_script="workflow/scripts/bssnp_bedGraph.pl",
        cpg="results/Illumina_pe/{protocol}/result_files/cpg.raw.vcf",
    output:
        "results/Illumina_pe/{protocol}/result_files/cpg.raw.CG.bedgraph",
    log:
        "logs/BisSNP_beGraph_{protocol}.log",
    conda:
        "../envs/bisSnp.yaml"
    shell:
        """
        perl  {input.perl_script} {input.cpg} CG
        """


rule rename_bissnp_output:
    input:
        "results/Illumina_pe/{protocol}/result_files/cpg.raw.CG.bedgraph",
    output:
        "results/Illumina_pe/{protocol}/result_files/bisSNP.bed",
    shell:
        "mv {input} {output}"
