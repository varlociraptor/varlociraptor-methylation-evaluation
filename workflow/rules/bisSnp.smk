# Newer versions of BisSNP (0.90.0 or 1.0.0) don't work
rule download_BisSNP:
    output:
        "../tools/Bis-tools/BisSNP-0.82.2.jar",
    log:
        "../logs/download_BisSNP.log",
    params:
        # base_dir="~/Documents/Promotion/",
        base_dir=config["base_dir"],
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p {params.base_dir}tools
        cd {params.base_dir}tools
        git clone https://github.com/dnaase/Bis-tools.git
        cd Bis-tools
        wget -O BisSNP-0.82.2.jar https://sourceforge.net/projects/bissnp/files/BisSNP-0.82.2/BisSNP-0.82.2.jar/download
        """


rule BisSNP:
    input:
        jar="../tools/Bis-tools/BisSNP-0.82.2.jar",
        genome="resources/genome.fasta",
        genome_index="resources/genome.fasta.fai",
        alignment="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/Illumina_pe/{protocol}/alignment_focused_downsampled_dedup_renamed.bam.bai",
    output:
        cpg="results/ref_tools/BisSNP/{protocol}/cpg.raw.vcf",
        snp="results/ref_tools/BisSNP/{protocol}/snp.raw.vcf",
        bed="results/ref_tools/BisSNP/{protocol}/cpg.raw.BisSNP-0.82.2.CG.bed",
    log:
        "logs/pb_CpG_tools_{protocol}.log",
    params:
        base_dir=config["base_dir"],
        prefix=lambda wildcards, input, output: os.path.splitext(output[0])[0].replace(
            ".combined", ""
        ),
        chromosome=chromosome_conf["chromosome"],
    shell:
        """
        java -Xmx10G -jar {input.jar} -R {input.genome} -T BisulfiteGenotyper -I {input.alignment} -vfn1 {output.cpg} -vfn2 {output.snp} -L {params.chromosome}
        perl ../tools/Bis-tools/utils/vcf2bed.pl {output.cpg} CG
        """
