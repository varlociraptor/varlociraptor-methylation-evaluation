# Fake data to simulate reads with Mason2


rule mason_download:
    output:
        mason_dir=directory("resources/tools/seqan/apps/mason2"),
        mason="resources/tools/seqan/apps/mason2/methylation_levels.h",
    log:
        "logs/mason/mason_download/download.log",
    conda:
        "../envs/shell_cmds.yaml"
    shell:
        """
        mkdir -p resources/tools 2> {log}
        cd resources/tools 2> {log}
        git clone git@github.com:seqan/seqan.git 2> {log}
        """


rule mason_fake_methylation:
    input:
        chrom="resources/{chrom}.fasta",
        index="resources/{chrom}.fasta.fai",
    output:
        methylation="resources/Simulate/simulated_data/{chrom}_meth.fa",
    log:
        "logs/mason/mason_fake_methylation/{chrom}.log",
    conda:
        "../envs/mason.yaml"
    params:
        seed=config["seed"],
    shell:
        """
        mkdir -p $(dirname {output.methylation})
        mason_methylation --in {input.chrom} \
            --methylation-levels \
            --meth-cg-sigma 0.3 \
            --meth-cg-mu 0.5 \
            --seed {params.seed} \
            --out {output.methylation}  2> {log}
        """


rule mason_fake_variants:
    input:
        chrom="resources/{chrom}.fasta",
    output:
        "resources/Simulate/simulated_data/{chrom}_variants.vcf",
    log:
        "logs/mason/mason_fake_variants/{chrom}.log",
    conda:
        "../envs/mason.yaml"
    params:
        seed=config["seed"],
    shell:
        """
        mason_variator --in-reference {input} \
            --out-vcf {output} \
            --seed {params.seed}  2> {log}
        """


rule mason_fake_reads:
    input:
        genome=expand(
            "resources/{chrom}.fasta",
            chrom=config["seq_platforms"].get("Simulate", []),
        ),
        genome_index=expand(
            "resources/{chrom}.fasta.fai",
            chrom=config["seq_platforms"].get("Simulate", []),
        ),
        variants=expand(
            "resources/Simulate/simulated_data/{chrom}_variants.vcf",
            chrom=config["seq_platforms"].get("Simulate", []),
        ),
        methylation=expand(
            "resources/Simulate/simulated_data/{chrom}_meth.fa",
            chrom=config["seq_platforms"].get("Simulate", []),
        ),
    output:
        f1="resources/Simulate/{sample}/{SRA}/{SRA}_1_trimmed.fastq.gz",
        f2="resources/Simulate/{sample}/{SRA}/{SRA}_2_trimmed.fastq.gz",
    log:
        "logs/mason/mason_fake_reads/{sample}_{SRA}.log",
    conda:
        "../envs/mason.yaml"
    threads: 16
    params:
        num_fragments=config.get("num_simulated_reads"),
    shell:
        """
        mason_simulator --input-reference {input.genome} \
                --input-vcf {input.variants} \
                --num-fragments {params.num_fragments} \
                --out {output.f1} \
                --out-right {output.f2} \
                --meth-fasta-in {input.methylation} \
                --enable-bs-seq \
                --num-threads {threads} \
                --illumina-read-length 150  2> {log}
        """


# Mason has a different meth ratio for forward and reverse strands.
# That is why we need to compute the coverage on the forward and reverse strand independently.
rule mason_alignment_forward:
    input:
        "resources/Simulate/simulated_data/alignment_focused_downsampled_dedup_renamed.bam",
    output:
        first="resources/Simulate/simulated_data/alignment_99.bam",
        second="resources/Simulate/simulated_data/alignment_147.bam",
        forward="resources/Simulate/simulated_data/alignment_forward.bam",
    log:
        "logs/mason/mason_alignment_forward/.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -b -f 64 -F 16 {input} > {output.first} 2> {log}
        samtools view -b -f 16 -F 64 {input} > {output.second} 2> {log}
        samtools merge {output.forward} {output.first} {output.second} 2> {log}
        """


rule mason_alignment_reverse:
    input:
        "resources/Simulate/simulated_data/alignment_focused_downsampled_dedup_renamed.bam",
    output:
        first="resources/Simulate/simulated_data/alignment_83.bam",
        second="resources/Simulate/simulated_data/alignment_163.bam",
        rev="resources/Simulate/simulated_data/alignment_reverse.bam",
    log:
        "logs/mason/mason_alignment_reverse/.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -b -f 16 -f 64 {input} > {output.first} 2> {log}
        samtools view -b -F 16 -F 64 {input} > {output.second} 2> {log}
        samtools merge {output.rev} {output.first} {output.second} 2> {log}
        """


rule mason_sort_oriented_reads:
    input:
        "resources/Simulate/simulated_data/alignment_{orientation}.bam",
    output:
        "resources/Simulate/simulated_data/alignment_sorted_{orientation}.bam",
    log:
        "logs/mason/mason_sort_oriented_reads/{orientation}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools sort -@ {threads}  {input} -o {output} 2> {log}"


rule mason_index_oriented_alignment:
    input:
        "resources/Simulate/simulated_data/alignment_sorted_{orientation}.bam",
    output:
        "resources/Simulate/simulated_data/alignment_sorted_{orientation}.bam.bai",
    log:
        "logs/mason/mason_index_oriented_alignment/{orientation}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input} 2> {log}"


# bcftools query -f '%CHROM\t%POS\t%REF\n' \
#   resources/J02459/candidates_1-of-1.bcf \
# | awk '{print $1 "\t" $2-1 "\t" $2-1+length($3)}' \
# > candidates.bed


rule candidates_to_bed:
    input:
        "resources/{chrom}/candidates.bcf",
    output:
        "resources/{chrom}/candidates.bed",
    log:
        "logs/mason/candidates_to_bed/{chrom}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%REF\n' {input} 2> {log} | \
        awk '{{print $1 "\t" $2-1 "\t" $2-1+length($3)}}' > {output}
        """


rule mason_coverage_orientation:
    input:
        bam="resources/Simulate/simulated_data/alignment_sorted_{orientation}.bam",
        bai="resources/Simulate/simulated_data/alignment_sorted_{orientation}.bam.bai",
        bed=expand(
            "resources/{chrom}/candidates.bed",
            chrom=(
                config["seq_platforms"].get("Simulate")
                if config["seq_platforms"].get("Simulate") != "genome"
                else "21"
            ),
        ),
    output:
        "resources/Simulate/simulated_data/{orientation}_cov.mosdepth.global.dist.txt",
        "resources/Simulate/simulated_data/{orientation}_cov.mosdepth.region.dist.txt",
        "resources/Simulate/simulated_data/{orientation}_cov.regions.bed.gz",
        summary="resources/Simulate/simulated_data/{orientation}_cov.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mason/mason_coverage_orientation/{orientation}.log",
    threads: 4  # This value - 1 will be sent to `--threads`
    params:
        extra="--no-per-base --use-median",  # optional
    wrapper:
        "v5.5.2/bio/mosdepth"


rule mason_candidates_vcf:
    input:
        "resources/{chrom}/candidates.bcf",
    output:
        "resources/{chrom}/candidates.vcf.gz",
    log:
        "logs/mason/mason_candidates_vcf/{chrom}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools view -o {output} {input} 2> {log}
        """


rule unzip:
    input:
        "resources/Simulate/simulated_data/{orientation}.regions.bed.gz",
    output:
        "resources/Simulate/simulated_data/{orientation}.regions.bed",
    log:
        "logs/mason/unzip/{orientation}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        gunzip -c {input} > {output}
        """


rule mason_compute_truth:
    input:
        cov_forward="resources/Simulate/simulated_data/forward_cov.regions.bed",
        cov_reverse="resources/Simulate/simulated_data/reverse_cov.regions.bed",
        methylation=(
            "resources/Simulate/simulated_data/"
            + (
                "genome"
                if config["seq_platforms"]["Simulate"] != "J02459"
                else "J02459"
            )
            + "_meth.fa"
        ),
        candidates="resources/{chrom}/candidates.vcf",
    output:
        "resources/Simulate/simulated_data/{chrom}_truth.csv",
    log:
        "logs/mason/mason_compute_truth/{chrom}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=32000,
    script:
        "../scripts/mason_ascii_to_meth.py"


rule mason_plot_truth_to_results:
    input:
        truth="resources/Simulate/simulated_data/{chrom}_truth.csv",
        results_rep="results/single_sample/Simulate/result_files/sample_df_simulated_data.parquet",
    output:
        report(
            "results/single_sample/Simulate/plots/simulated_data_{chrom}.{plot_type}",
            category="single_sample",
            subcategory="Simulated",
            labels=lambda wildcards: {
                "file_type": "heatmap",
                "sample": f"simulated_data",
            },
            caption="../report/heatmap.rst",
        ),
    log:
        "logs/mason/mason_plot_truth_to_results/{chrom}_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get("Simulate", [])
        + [f"varlo_{fdr}" for fdr in config["fdr_alpha"]],
        bin_size=lambda wildcards: config["heatmap_bin_size"],
    script:
        "../scripts/plot_mason_results.py"


rule compute_precision_recall:
    input:
        truth="resources/Simulate/simulated_data/{chrom}_truth.csv",
        results_rep="results/single_sample/Simulate/result_files/sample_df_simulated_data.parquet",
        no_bias="results/single_sample/Simulate/called/simulated_data_no_bias/result_files/varlo_0.01.parquet",
        # coverage="resources/Simulate/simulated_data/complete_cov.regions.bed",
        # tool="results/{platform}/{protocol}/result_files/{method}.parquet",
    output:
        precall="results/single_sample/Simulate/plots/precall_{chrom}.{plot_type}",
        cov_dist="results/single_sample/Simulate/plots/{chrom}_cov_dist.{plot_type}",
    log:
        "logs/mason/compute_precision_recall/{chrom}_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    params:
        meth_callers=lambda wildcards: config["ref_tools"].get("Simulate", [])
        + [f"varlo_{fdr}" for fdr in config["fdr_alpha"]],
    script:
        "../scripts/compute_precision_recall.py"


# Truth positions whose coverage falls into the given coverage bin
rule truth_positions_per_coverage_bin:
    input:
        "resources/Simulate/simulated_data/{chrom}_truth.csv",
    output:
        "resources/Simulate/simulated_data/{chrom}_coverage_bin_{bin}.tsv",
    log:
        "logs/mason/truth_positions_per_coverage_bin/{chrom}_{bin}.log",
    wildcard_constraints:
        bin=r"\d+",
    conda:
        "../envs/python.yaml"
    params:
        coverage_bins=config["calibration"]["coverage_bins"],
    script:
        "../scripts/truth_positions_per_coverage_bin.py"


# Split the simulated varlociraptor calls into coverage bins defined by the truth coverage.
# The scattered calls are merged first, since global FDR control has to see all calls of a bin at once.
rule split_calls_by_coverage_bin:
    input:
        calls=gather.split_candidates(
            "results/single_sample/Simulate/fdr_calibration/simulated_data/calls_{scatteritem}.bcf"
        ),
        positions="resources/Simulate/simulated_data/{chrom}_coverage_bin_{bin}.tsv",
    output:
        "results/single_sample/Simulate/fdr_calibration/simulated_data/coverage_bins/{chrom}/bin_{bin}.bcf",
    log:
        "logs/mason/split_calls_by_coverage_bin/{chrom}_{bin}.log",
    wildcard_constraints:
        bin=r"\d+",
    conda:
        "../envs/samtools.yaml"
    shell:
        "(bcftools concat -Ou {input.calls} | bcftools view -T {input.positions} -Ob -o {output}) 2> {log}"


# Select the calls of a coverage bin with varlociraptor's global FDR control on the HIGH event
rule control_fdr_per_coverage_bin:
    input:
        varlo=varlociraptor_binary(),
        calls="results/single_sample/Simulate/fdr_calibration/simulated_data/coverage_bins/{chrom}/bin_{bin}.bcf",
    output:
        "results/single_sample/Simulate/fdr_calibration/simulated_data/coverage_bins/{chrom}/bin_{bin}_fdr_{alpha}.bcf",
    log:
        "logs/mason/control_fdr_per_coverage_bin/{chrom}_{bin}_{alpha}.log",
    wildcard_constraints:
        bin=r"\d+",
    params:
        varlo=varlociraptor_cmd,
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{params.varlo} filter-calls control-fdr {input.calls} --events HIGH --fdr {wildcards.alpha} --mode global-strict > {output} 2> {log}"


rule plot_fdr_calibration:
    input:
        truth="resources/Simulate/simulated_data/{chrom}_truth.csv",
        # Calls per coverage bin before FDR control, to count the sites per bin
        unfiltered=expand(
            "results/single_sample/Simulate/fdr_calibration/simulated_data/coverage_bins/{{chrom}}/bin_{bin}.vcf",
            bin=range(len(config["calibration"]["coverage_bins"])),
        ),
        filtered=expand(
            "results/single_sample/Simulate/fdr_calibration/simulated_data/coverage_bins/{{chrom}}/bin_{bin}_fdr_{alpha}.vcf",
            bin=range(len(config["calibration"]["coverage_bins"])),
            alpha=config["calibration"]["fdr_levels"],
        ),
    output:
        plot=report(
            "results/single_sample/Simulate/plots/fdr_calibration_{chrom}.{plot_type}",
            category="single_sample",
            subcategory="Simulated",
            labels={
                "file_type": "fdr_calibration",
                "sample": "simulated_data",
            },
            caption="../report/fdr_calibration.rst",
        ),
        table="results/single_sample/Simulate/plots/fdr_calibration_{chrom}_{plot_type}.csv",
    log:
        "logs/mason/plot_fdr_calibration/{chrom}_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    params:
        fdr_levels=config["calibration"]["fdr_levels"],
        # Coverage stratifications
        coverage_bins=config["calibration"]["coverage_bins"],
        # Methylation rate in the truth above which a site is defined as HIGH
        high_threshold=50,
    script:
        "../scripts/plot_fdr_calibration.py"


rule plot_pit_calibration:
    input:
        truth="resources/Simulate/simulated_data/{chrom}_truth.csv",
        calls="results/single_sample/Simulate/fdr_calibration/simulated_data/result_files/varlo.bed",
    output:
        plot=report(
            "results/single_sample/Simulate/plots/pit_calibration_{chrom}.{plot_type}",
            category="single_sample",
            subcategory="Simulated",
            labels={
                "file_type": "pit_calibration",
                "sample": "simulated_data",
            },
            caption="../report/pit_calibration.rst",
        ),
        table="results/single_sample/Simulate/plots/pit_calibration_{chrom}_{plot_type}.csv",
    log:
        "logs/mason/plot_pit_calibration/{chrom}_{plot_type}.log",
    conda:
        "../envs/python.yaml"
    params:
        # Coverage stratifications
        coverage_bins=config["calibration"]["coverage_bins"],
        num_bins=25,
        # Seed for the randomized PIT at the point mass θ = 0
        seed=config["seed"],
    script:
        "../scripts/plot_pit_calibration.py"
