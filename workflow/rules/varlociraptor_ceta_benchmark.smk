
rule download_varlociraptor:
    output:
        direc=directory(
            "resources/tools/ceta_comp/varlociraptor",
        ),
        execu="resources/tools/ceta_comp/varlociraptor/target/debug/varlociraptor",
    log:
        "logs/download_varlociraptor.log",
    shell:
        """
        PARENT_DIR=$(dirname {output.direc})        
        mkdir -p $PARENT_DIR        
        cd $PARENT_DIR
        git clone https://github.com/varlociraptor/varlociraptor.git
        cd varlociraptor
        git checkout ceta_benchmarking
        cargo build --release
        """


# candidates with variant:    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001542345.1/
rule download_untreated_fastqs:
    output:
        fq1="resources/Illumina_pe/untreated/dummy/dummy_1.fastq",
        fq2="resources/Illumina_pe/untreated/dummy/dummy_2.fastq",
    params:
        fq1_url="https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R1.fastq.gz",
        fq2_url="https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R2.fastq.gz",
    shell:
        """
        mkdir -p $(dirname {output.fq1})

        # Download + Entpacken Read 1
        curl -L -C - {params.fq1_url} | gunzip -c > {output.fq1}.tmp && mv {output.fq1}.tmp {output.fq1}

        # Download + Entpacken Read 2
        curl -L -C - {params.fq2_url} | gunzip -c > {output.fq2}.tmp && mv {output.fq2}.tmp {output.fq2}
        """


rule trim_untreated_fastqs:
    input:
        fq1="resources/Illumina_pe/untreated/dummy/dummy_1.fastq",
        fq2="resources/Illumina_pe/untreated/dummy/dummy_2.fastq",
    output:
        fq1="resources/Illumina_pe/untreated/dummy/dummy_1_trimmed.fastq",
        fq2="resources/Illumina_pe/untreated/dummy/dummy_2_trimmed.fastq",
    log:
        "logs/data/trim_untreated_fastqs.log",
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp --in1 {input.fq1} --in2 {input.fq2} --out1 {output.fq1} --out2 {output.fq2} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis 2> {log}"


# TODO: This is ambiguous with rule bwameth_index
rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        idx=multiext("resources/genome.{alg}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/resources/genome.{alg}.log",
    params:
        extra=lambda w: f"-a {w.alg}",
    wildcard_constraints:
        alg="(?!bwameth)",
    wrapper:
        "v5.10.0/bio/bwa/index"


rule bwa_mem:
    input:
        reads=[
            "resources/Illumina_pe/untreated/dummy/dummy_1_trimmed.fastq",
            "resources/Illumina_pe/untreated/dummy/dummy_2_trimmed.fastq",
        ],
        idx=multiext("resources/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "resources/Illumina_pe/untreated/dummy/alignment.bam",
    log:
        "logs/bwa_mem.log",
    params:
        sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v7.6.0/bio/bwa/mem"


rule download_variant_truth:
    output:
        "resources/ceta/candidates_variants.vcf.gz",
    params:
        url="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    log:
        "logs/variants/download_variant_truth.log",
    shell:
        "wget -q -O {output} {params.url} 2> {log}"


rule index_variant_truth:
    input:
        "resources/ceta/candidates_variants.vcf.gz",
    output:
        "resources/ceta/candidates_variants.vcf.gz.csi",
    log:
        "logs/variants/index_variant_truth.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools index -f {input} 2> {log}"


rule variants_focus_chromosome:
    input:
        vcf="resources/ceta/candidates_variants.vcf.gz",
        idx="resources/ceta/candidates_variants.vcf.gz.csi",
    output:
        "resources/ceta/variants_focus_chromosome.bcf",
    threads: 4
    shell:
        """
        bcftools view --regions chr21 {input.vcf} \
            --threads {threads} \
            -Ob -o {output}
        bcftools index -f {output}
        """


rule rename_variant_chromosome:
    input:
        "resources/ceta/variants_focus_chromosome.bcf",
    output:
        "resources/ceta/candidates_variants_renamed.bcf",
    threads: 4
    log:
        "logs/variants/rename_variant_chromosome.log",
    shell:
        """
        bcftools annotate --rename-chrs <(echo -e "chr21\t21") -o {output} -O b {input} 2> {log}
        """


rule split_ceta_candidates:
    input:
        "resources/ceta/candidates_variants_renamed.bcf",
    output:
        scatter.split_candidates("resources/ceta/candidates_{scatteritem}.bcf"),
    log:
        "logs/candidates/split_ceta_candidates.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output} 2> {log}"


rule varlociraptor_ceta_preprocess:
    input:
        varlo="resources/tools/ceta_comp/varlociraptor/target/debug/varlociraptor",
        chromosome="resources/chromosome_21.fasta",
        genome_index="resources/chromosome_21.fasta.fai",
        alignments="resources/Illumina_pe/{sample}/candidate_specific/alignment_{scatteritem}.bam",
        alignment_index="resources/Illumina_pe/{sample}/candidate_specific/alignment_{scatteritem}.bam.bai",
        candidates="resources/ceta/candidates_{scatteritem}.bcf",
    output:
        "results/ceta_benchmark/preprocessed/Illumina_pe/{sample}/normal_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/Illumina_pe/{sample}/compute_meth_observations_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    resources:
        mem_mb=16000,
    shell:
        """
        {input.varlo} preprocess variants --omit-mapq-adjustment {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --max-depth 5000 > {output} 2> {log}
        """


rule varlociraptor_ceta_call_single:
    input:
        varlo="resources/tools/ceta_comp/varlociraptor/target/debug/varlociraptor",
        preprocess_obs="results/ceta_benchmark/preprocessed/Illumina_pe/{sample}/normal_{scatteritem}.bcf",
        scenario=lambda wc: (
            f"resources/scenarios/ceta_benchmarks/scenario_untreated.yaml"
            if wc.sample == "untreated"
            else f"resources/scenarios/ceta_benchmarks/scenario_converted.yaml"
        ),
    output:
        "results/ceta_benchmark/Illumina_pe/called/{sample}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/Illumina_pe/called/{sample}/call_methylation_{scatteritem}.log",
    wildcard_constraints:
        sample="(?!ceta_multi).*",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output} 2> {log}"


rule varlociraptor_ceta_call_multi_emseq:
    input:
        varlo="resources/tools/ceta_comp/varlociraptor/target/debug/varlociraptor",
        emseq="results/ceta_benchmark/preprocessed/Illumina_pe/EMSeq_HG002_LAB01_REP01/normal_{scatteritem}.bcf",
        untreated="results/ceta_benchmark/preprocessed/Illumina_pe/untreated/normal_{scatteritem}.bcf",
        scenario="resources/scenarios/ceta_benchmarks/scenario_common.yaml",
    output:
        "results/ceta_benchmark/Illumina_pe/called/ceta_multi/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/multi_sample/untreated_emseq/call_methylation_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs  emseq={input.emseq}  untreated={input.untreated} > {output} 2> {log}"

rule varlociraptor_ceta_call_multi_all:
    input:
        varlo="resources/tools/ceta_comp/varlociraptor/target/debug/varlociraptor",
        emseq="results/ceta_benchmark/preprocessed/Illumina_pe/EMSeq_HG002_LAB01_REP01/normal_{scatteritem}.bcf",
        untreated="results/ceta_benchmark/preprocessed/Illumina_pe/untreated/normal_{scatteritem}.bcf",
        methylseq="results/ceta_benchmark/preprocessed/Illumina_pe/MethylSeq_HG002_LAB01_REP01/normal_{scatteritem}.bcf",
        scenario="resources/scenarios/ceta_benchmarks/scenario_common_all.yaml",
    output:
        "results/ceta_benchmark/Illumina_pe/called/ceta_multi_all/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/multi_sample/untreated_emseq/call_methylation_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs  emseq={input.emseq} methylseq={input.methylseq} untreated={input.untreated} > {output} 2> {log}"


rule event_probs_df:
    input:
        tool="results/ceta_benchmark/Illumina_pe/called/{sample}/result_files/varlo.bed",
        cg_candidates="resources/21/candidates.bcf",
    output:
        "results/ceta_benchmark/Illumina_pe/called/{sample}/result_files/events_{fdr}.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plot_results/event_probs_df/{sample}_{fdr}.log",
    params:
        alpha=lambda wildcards: wildcards.fdr,
    resources:
        mem_mb=64000,
    script:
        "../scripts/event_probs_df.py"


rule plot_ceta_probs:
    input:
        "results/ceta_benchmark/Illumina_pe/called/MethylSeq_HG002_LAB01_REP01/result_files/events_{fdr}.parquet",
        "results/ceta_benchmark/Illumina_pe/called/EMSeq_HG002_LAB01_REP01/result_files/events_{fdr}.parquet",
        "results/ceta_benchmark/Illumina_pe/called/ceta_multi/result_files/events_{fdr}.parquet",
        "results/ceta_benchmark/Illumina_pe/called/ceta_multi_all/result_files/events_{fdr}.parquet",
        "results/ceta_benchmark/Illumina_pe/called/untreated/result_files/events_{fdr}.parquet",
    output:
        "results/ceta_benchmark/Illumina_pe/called/result_files/combined_{fdr}.html",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/plot_ceta_probs_{fdr}.log",
    resources:
        mem_mb=16000,
    script:
        "../scripts/plot_ceta_probs.py"
