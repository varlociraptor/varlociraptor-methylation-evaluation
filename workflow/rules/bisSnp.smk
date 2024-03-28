# Newer versions of BisSNP (0.90.0 or 1.0.0) don't work
rule download_BisSNP:
    output:
        "resources/ref_tools/Bis-tools/BisSNP-0.82.2.jar",
    log:
        "../logs/download_BisSNP.log",
    params:
        # base_dir="~/Documents/Promotion/",
        pipeline_path=config["pipeline_path"],
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p {params.pipeline_path}resources/ref_tools
        cd {params.pipeline_path}/resources/ref_tools
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
        cpg="results/ref_tools/BisSNP/{protocol}/cpg.raw.vcf",
        snp="results/ref_tools/BisSNP/{protocol}/snp.raw.vcf",
    log:
        "logs/pb_CpG_tools_{protocol}.log",
    conda:
        "../envs/bisSnp.yaml"
    params:
        base_dir=config["base_dir"],
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ""
        ),
        chromosome=chromosome_conf["chromosome"],
    shell:
        """
        java -Xmx10G -jar {input.jar} -R {input.genome} -T BisulfiteGenotyper -I {input.alignment} -vfn1 {output.cpg} -vfn2 {output.snp} -L {params.chromosome}
        """


# We do not use the official perl script in resources/ref_tools/Bis-tools/utils/vcf2bedGraph.pl because it does not work
# We copied the script and removed line 79: next unless ($splitin[6] eq "PASS" || $splitin[6] eq "Infinity"); because it never triggers
rule BisSNP_bedGraph:
    input:
        perl_script="workflow/scripts/bssnp_bedGraph.pl",
        cpg="results/ref_tools/BisSNP/{protocol}/cpg.raw.vcf",
    output:
        bed="results/ref_tools/BisSNP/{protocol}/cpg.raw.CG.bedgraph",
    log:
        "logs/BisSNP_beGraph_{protocol}.log",
    conda:
        "../envs/bisSnp.yaml"
    shell:
        """
        perl  {input.perl_script} {input.cpg} CG
        """
