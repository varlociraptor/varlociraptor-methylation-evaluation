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
        genome=lambda wildcards: (
            expand(
                "resources/{chrom}.fasta",
                chrom=config["seq_platforms"].get(wildcards.platform),
            )
            if wildcards.sample.startswith("simulated_data")
            else ["resources/genome.fasta"]
        ),
        genome_index=lambda wildcards: (
            expand(
                "resources/{chrom}.fasta.fai",
                chrom=config["seq_platforms"].get(wildcards.platform),
            )
            if wildcards.sample.startswith("simulated_data")
            else ["resources/genome.fasta.fai"]
        ),
        # This maybe for simulated  data:
        alignment="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
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
        # Regenerate the FASTA index to ensure it matches the copied genome
        samtools faidx {output.genome} 2>> {log}
        samtools view -H {input.alignment} > /tmp/header.sam 2>> {log}
        # Add read group header if it doesn't exist
        if ! grep -q "^@RG" /tmp/header.sam; then
            echo "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}" >> /tmp/header.sam
        fi
        samtools reheader /tmp/header.sam {input.alignment} | samtools sort -o {output.alignment} - 2>> {log}
        samtools index {output.alignment} 2>> {log}
        rm /tmp/header.sam
        """



rule bissnp_extract:
    input:
        jar="resources/ref_tools/Bis-tools/{platform}/{sample}/BisSNP-0.82.2.jar",
        genome="resources/ref_tools/Bis-tools/{platform}/{sample}/genome.fasta",
        # genome_dict="resources/ref_tools/Bis-tools/{platform}/{sample}/genome.dict",
        genome_index="resources/ref_tools/Bis-tools/{platform}/{sample}/genome.fasta.fai",
        alignment="resources/ref_tools/Bis-tools/{platform}/{sample}/alignment.bam",
        alignment_index="resources/ref_tools/Bis-tools/{platform}/{sample}/alignment.bam.bai",
    output:
        cpg="results/single_sample/{platform}/called/{sample}/result_files/bissnp_cpg.raw.vcf",
        snp="results/single_sample/{platform}/called/{sample}/result_files/bissnp_snp.raw.vcf",
    conda:
        "../envs/openjdk.yaml"
    params:
        loc_flag=lambda wildcards: "" if chromosome_by_seq_platform.get(wildcards.platform) == "genome" else f"-L {chromosome_by_seq_platform.get(wildcards.platform)}",
    log:
        "logs/bissnp/bissnp_extract/{platform}_{sample}.log",
    benchmark:
        repeat("benchmarks/{platform}/bisSNP/bissnp_extract/{sample}.txt", config["benchmark_repeats"])
    threads: 8
    resources:
        mem_mb=16000,
    shell:
        """
        stdbuf -oL -eL java -Xmx10G -jar {input.jar} \
            -R {input.genome} \
            -nt {threads} \
            -T BisulfiteGenotyper \
            -I {input.alignment} \
            -vfn1 {output.cpg} \
            -vfn2 {output.snp} \
            {params.loc_flag} > {log} 2>&1
        """

# We do not use the official perl script in resources/ref_tools/Bis-tools/utils/vcf2bedGraph.pl because it does not work
# We copied the script and removed line 79: next unless ($splitin[6] eq "PASS" || $splitin[6] eq "Infinity"); because it never triggers
# Furthermore we added coverage information to the output:
#   - line 68: my $head_line = "track type=bedGraph name=${cpg_name_output}.${bissnp_version} description=\"$type methylation level and coverage\" visibility=3";
#   - line 104: my $out_line = "$chr\t$start\t$end\t$methy\t$ct_reads";
rule bissnp_create_bedgraph:
    input:
        perl_script=workflow.source_path("../scripts/bissnp_bedGraph.pl"),
        cpg="results/single_sample/{platform}/called/{sample}/result_files/bissnp_cpg.raw.vcf",
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
        candidates=lambda wildcards: expand(
            "resources/{chrom}/candidates.bcf",
            chrom=config["seq_platforms"].get(wildcards.platform),
        ),
        candidates_index=lambda wildcards: expand(
            "resources/{chrom}/candidates.bcf.csi",
            chrom=config["seq_platforms"].get(wildcards.platform),
        ),
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/bisSNP.bed",
    log:
        "logs/bissnp/bissnp_merge_positions/{platform}_{sample}.log",
    conda:
        "../envs/pysam.yaml"
    resources:
        mem_mb=64000
    script:
        "../scripts/merge_forward_reverse_positions.py"
