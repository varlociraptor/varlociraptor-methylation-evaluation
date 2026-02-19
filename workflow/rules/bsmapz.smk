# Clone BSMAPz and compile bsmapz locally
# The official conda environment does not work properly (Illegal instruction     (core dumped))
rule bsmapz_clone_and_build:
    output:
        binary="resources/ref_tools/BSMAPz/bsmapz",
        meth_extractor="resources/ref_tools/BSMAPz/methratio.py",
    log:
        "logs/bsmapz/bsmapz_clone_and_build/download.log",
    conda:
        "../envs/general.yaml"
    shell:
        """
        mkdir -p resources/ref_tools
        cd resources/ref_tools
        if [ ! -d BSMAPz ]; then
        git clone https://github.com/zyndagj/BSMAPz.git
        fi
        cd BSMAPz
        make bsmapz
        echo "BSMAPz log output"
        ls
        pwd
        cd ../..
        cp bsmapz ../../bsmapz
        """


# # Download the newer methylation extractor from BSMAPz
# rule bsmap_download_methratio:
#     output:
#         "resources/ref_tools/BSMAPz/methratio.py",
#     log:
#         "logs/bsmapz/downlzoad_methratio.log",
#     shell:
#         """
#         mkdir -p $(dirname {output})
#         wget -O {output} https://raw.githubusercontent.com/zyndagj/BSMAPz/master/methratio.py > {log} 2>&1
# """


# Run BSMAPz to compute methylation alignments
rule bsmapz_compute_meth:
    input:
        genome=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"].get(wildcards.platform),
        ),
        alignment="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        bsmapz_binary="resources/ref_tools/BSMAPz/bsmapz",
    output:
        temp("results/single_sample/{platform}/called/{sample}/result_files/out.sam"),
    log:
        "logs/bsmapz/bsmapz_compute/{platform}_{sample}.log",
    resources:
        mem_mb=32000,
    benchmark:
        "benchmarks/{platform}/bsmap/bsmap_compute/{sample}.bwa.benchmark.txt"
    conda:
        "../envs/general.yaml"
    threads: 1
    shell:
        """

        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output})
        chmod +x {input.bsmapz_binary}
        {input.bsmapz_binary} -a {input.alignment} -d {input.genome} -o {output} -p {threads} -w 100 -v 0.07 -m 50 -x 300 > {log} 2>&1
        """


# Extract methylation ratios
rule bsmapz_extract:
    input:
        genome=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["seq_platforms"].get(wildcards.platform),
        ),
        genome_index=lambda wildcards: expand(
            "resources/chromosome_{chrom}.fasta.fai",
            chrom=config["seq_platforms"].get(wildcards.platform),
        ),
        bsmap_sam="results/single_sample/{platform}/called/{sample}/result_files/out.sam",
        meth_extractor="resources/ref_tools/BSMAPz/methratio.py",
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/methylation_ratios.bed",
    log:
        "logs/bsmapz/bsmapz_extract/{platform}_{sample}.log",
    params:
        chromosome=lambda wildcards: chromosome_by_seq_platform.get(wildcards.platform, "21"),
    conda:
        "../envs/bsmapz.yaml"
    benchmark:
        "benchmarks/{platform}/bsmap/bsmap_extract/{sample}.bwa.benchmark.txt"
    shell:
        "python {input.meth_extractor} -c={params.chromosome} --ref={input.genome} --out={output} {input.bsmap_sam} -g -x CG 2> {log}"


# Rename output file to standardized name
rule bsmapz_rename_output:
    input:
        "results/single_sample/{platform}/called/{sample}/result_files/methylation_ratios.bed",
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/bsMap.bed",
    log:
        "logs/bsmapz/bsmapz_rename_output/{platform}_{sample}.log",
    conda:
        "../envs/general.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        cp {input} {output} 2> {log}
        rm {input} 2> {log}
        """
