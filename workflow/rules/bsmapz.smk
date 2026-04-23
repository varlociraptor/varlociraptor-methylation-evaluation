# Clone BSMAPz and compile bsmapz locally
# The official conda environment does not work properly (Illegal instruction     (core dumped))
# rule bsmapz_clone_and_build:
#     output:
#         binary="resources/ref_tools/BSMAPz/bsmapz",
#         meth_extractor="resources/ref_tools/BSMAPz/methratio.py",
#     log:
#         "logs/bsmapz/bsmapz_clone_and_build/download_bsmapz.log",
#     conda:
#         "../envs/general.yaml"
#     shell:
#         """
#         build_dir=$(mktemp -d)

#         git clone https://github.com/zyndagj/BSMAPz.git $build_dir/BSMAPz 2> {log}

#         export CFLAGS="-I$CONDA_PREFIX/include"
#         export CXXFLAGS="-I$CONDA_PREFIX/include"
#         export LDFLAGS="-L$CONDA_PREFIX/lib"

#         make -C "$build_dir/BSMAPz" bsmapz >> {log} 2>&1

#         cp "$build_dir/BSMAPz/bsmapz"      {output.binary}
#         cp "$build_dir/BSMAPz/methratio.py" {output.meth_extractor}

#         rm -rf "$build_dir"
#         """



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
            "resources/{chrom}.fasta",
            chrom=config["seq_platforms"].get(wildcards.platform),
        ),
        alignment="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam",
        alignment_index="resources/{platform}/{sample}/alignment_focused_downsampled_dedup_renamed.bam.bai",
        # bsmapz_binary="resources/ref_tools/BSMAPz/bsmapz",
    output:
        temp("results/single_sample/{platform}/called/{sample}/result_files/out.unsorted.bam"),
    log:
        "logs/bsmapz/bsmapz_compute/{platform}_{sample}.log",
    resources:
        mem_mb=16000,
    benchmark:
        repeat("benchmarks/{platform}/bsmap/bsmap_compute/{sample}.bwa.benchmark.txt", config["benchmark_repeats"])
    conda:
        "../envs/bsmapz.yaml"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output})
        bsmapz -a {input.alignment} -d {input.genome} -o {output} -p {threads} -w 100 -v 0.07 -m 50 -x 300 > {log} 2>&1
        """

# Sort BSMAPz output BAM by coordinate (BSMAPz does not guarantee sorted output)
rule bsmapz_sort_out_bam:
    input:
        "results/single_sample/{platform}/called/{sample}/result_files/out.unsorted.bam",
    output:
        temp("results/single_sample/{platform}/called/{sample}/result_files/out.bam"),
    log:
        "logs/bsmapz/bsmapz_sort_out_bam/{platform}_{sample}.log",
    resources:
        mem_mb=8000,
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2> {log}"



# Index out.bam for region-based splitting
rule bsmapz_index_out_bam:
    input:
        "results/single_sample/{platform}/called/{sample}/result_files/out.bam",
    output:
        temp(
            "results/single_sample/{platform}/called/{sample}/result_files/out.bam.bai"
        ),
    log:
        "logs/bsmapz/bsmapz_index_out_bam/{platform}_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input} 2> {log}"


# Split out.bam by candidate region for parallel methylation extraction
rule bsmapz_extract_scatter_bam:
    input:
        alignment="results/single_sample/{platform}/called/{sample}/result_files/out.bam",
        index="results/single_sample/{platform}/called/{sample}/result_files/out.bam.bai",
        candidate=lambda wildcards: f"resources/{chromosome_by_seq_platform.get(wildcards.platform)}/candidates_{wildcards.scatteritem}.bed",
    output:
        temp(
            "results/single_sample/{platform}/called/{sample}/result_files/out_{scatteritem}.bam"
        ),
    log:
        "logs/bsmapz/bsmapz_extract_scatter_bam/{platform}_{sample}_{scatteritem}.log",
    params:
        chromosome=lambda wildcards: chromosome_by_seq_platform.get(
            wildcards.platform, "21"
        ),
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -b -L {input.candidate} {input.alignment} > {output} 2> {log}

        if [ $(samtools view -c {output}) -eq 0 ]; then
            samtools view -H {input.alignment} > {output}.temp.sam
            samtools view {input.alignment} | tail -n 1 >> {output}.temp.sam
            samtools view -bS {output}.temp.sam > {output}
            rm {output}.temp.sam
        fi
        """


# Extract methylation ratios per region (scattered)
rule bsmapz_extract:
    input:
        genome=lambda wildcards: expand(
            "resources/{chrom}.fasta",
            chrom=config["seq_platforms"].get(wildcards.platform),
        ),
        genome_index=lambda wildcards: expand(
            "resources/{chrom}.fasta.fai",
            chrom=config["seq_platforms"].get(wildcards.platform),
        ),
        bsmap_bam="results/single_sample/{platform}/called/{sample}/result_files/out_{scatteritem}.bam",
        meth_extractor="resources/ref_tools/BSMAPz/methratio.py",
    output:
        temp(
            "results/single_sample/{platform}/called/{sample}/result_files/methylation_ratios_{scatteritem}.bed"
        ),
    log:
        "logs/bsmapz/bsmapz_extract/{platform}_{sample}_{scatteritem}.log",
    params:
        chromosome_flag=lambda wildcards: f"-c={chromosome_by_seq_platform.get(wildcards.platform)}"
            if chromosome_by_seq_platform.get(wildcards.platform) != "genome"
                else ""
    conda:
        "../envs/bsmapz.yaml"
    resources:
        mem_mb=64000
    benchmark:
        repeat("benchmarks/{platform}/bsmap/bsmap_extract/{sample}_{scatteritem}.bwa.benchmark.txt", config["benchmark_repeats"])
    shell:
        "python {input.meth_extractor} {params.chromosome_flag} --ref={input.genome[0]} --out={output} {input.bsmap_bam} -g -x CG 2> {log}"


# Gather scattered methylation ratio BEDs into one file
rule bsmapz_gather_methylation:
    input:
        gather.split_candidates(
            "results/single_sample/{{platform}}/called/{{sample}}/result_files/methylation_ratios_{scatteritem}.bed"
        ),
    output:
        "results/single_sample/{platform}/called/{sample}/result_files/methylation_ratios.bed",
    log:
        "logs/bsmapz/bsmapz_gather_methylation/{platform}_{sample}.log",
    conda:
        "../envs/general.yaml"
    shell:
        """
        head -n1 $(echo {input} | tr ' ' '\n' | head -n1) > {output} 2> {log}
        for f in {input}; do tail -n +2 "$f"; done >> {output} 2>> {log}
        """


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
