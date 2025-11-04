# candidates with variant:    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/

rule download_untreated_fastqs:
    output:
        fq1="resources/Illumina_pe/UNTREATED/dummy/dummy_1.fastq",
        fq2="resources/Illumina_pe/UNTREATED/dummy/dummy_2.fastq",
    params:
        fq1_url="https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R1.fastq.gz",
        fq2_url="https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R2.fastq.gz"
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
        fq1="resources/Illumina_pe/UNTREATED/dummy/dummy_1.fastq",
        fq2="resources/Illumina_pe/UNTREATED/dummy/dummy_2.fastq",
    output:
        fq1="resources/Illumina_pe/UNTREATED/dummy/dummy_1_trimmed.fastq",
        fq2="resources/Illumina_pe/UNTREATED/dummy/dummy_2_trimmed.fastq",
    log:
        "logs/data/trim_untreated_fastqs.log",
    conda:
        "../envs/fastp.yaml"
    wildcard_constraints:
        sample="^(?!simulated_data$).*",
    shell:
        "fastp --in1 {input.fq1} --in2 {input.fq2} --out1 {output.fq1} --out2 {output.fq2} --length_required 2 --disable_quality_filtering -z 4 --trim_poly_g --overrepresentation_analysis 2> {log}"

# TODO: This is ambiguous with rule bwameth_index
rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        idx=multiext("resources/{genome}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/{genome}.log",
    wrapper:
        "v5.10.0/bio/bwa/index"

rule bwa_mem:
    input:
        reads=["resources/Illumina_pe/UNTREATED/dummy/dummy_1_trimmed.fastq", "resources/Illumina_pe/UNTREATED/dummy/dummy_2_trimmed.fastq"],
        idx=multiext("resources/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "resources/Illumina_pe/UNTREATED/dummy/alignment.bam",
    log:
        "logs/bwa_mem.log",
    params:
        sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v7.6.0/bio/bwa/mem"





rule download_variant_bcf:
    output:
        bcf="resources/ceta/candidates_variants.bcf"
    params:
        url="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    threads: 4
    shell:
        """
        # temporäre Dateien
        tmp_vcf=$(mktemp --suffix=.vcf.gz)

        # Datei herunterladen
        wget -q -O $tmp_vcf {params.url}

        # chr21 extrahieren und als BCF speichern
        bcftools view --regions chr21 $tmp_vcf \
            --threads {threads} \
            -Ob -o {output.bcf}

        # Index erzeugen (optional, aber sehr sinnvoll)
        bcftools index -f {output.bcf}
        """



rule split_ceta_candidates:
    input:
        "resources/ceta/candidates_variants.bcf",
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
        chromosome=            "resources/chromosome_21.fasta",
        genome_index=            "resources/chromosome_21.fasta.fai",
        alignments="resources/Illumina_pe/{sample}/candidate_specific/alignment_{scatteritem}.bam",
        alignment_index="resources/Illumina_pe/{sample}/candidate_specific/alignment_{scatteritem}.bam.bai",
        candidates=            "resources/ceta/candidates_{scatteritem}.bcf",
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
        {input.varlo} preprocess variants --omit-mapq-adjustment {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --read-type Illumina_pe --max-depth 5000 > {output} 2> {log}
        """

rule varlociraptor_ceta_call_single:
    input:
        varlo="resources/tools/ceta_comp/varlociraptor/target/debug/varlociraptor",
        preprocess_obs="results/ceta_benchmark/preprocessed/Illumina_pe/{sample}/normal_{scatteritem}.bcf",
        scenario="resources/scenarios/scenario.yaml",
    output:
        "results/ceta_benchmark/Illumina_pe/called/{sample}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/Illumina_pe/called/{sample}/call_methylation_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    wildcard_constraints:
        seq_platform="(?!multi_sample).*",
        sample="(?!ceta_multi).*",
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output} 2> {log}"



rule varlociraptor_ceta_call_multi:
    input:
        varlo="resources/tools/ceta_comp/varlociraptor/target/debug/varlociraptor",
        emseq="results/ceta_benchmark/preprocessed/Illumina_pe/EMSeq_HG002_LAB01_REP01/normal_{scatteritem}.bcf",
        untreat="results/ceta_benchmark/preprocessed/Illumina_pe/UNTREATED/normal_{scatteritem}.bcf",
        scenario="resources/scenarios/scenario_ceta_benchmark.yaml",
    output:
        "results/ceta_benchmark/Illumina_pe/called/ceta_multi/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/multi_sample/untreat_emseq/call_methylation_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs  emseq={input.emseq}  untreat={input.untreat} > {output} 2> {log}"





rule plot_ceta_probs:
    input:
        "results/ceta_benchmark/Illumina_pe/called/MethylSeq_HG002_LAB01_REP01/result_files/varlo.parquet",
        "results/ceta_benchmark/Illumina_pe/called/EMSeq_HG002_LAB01_REP01/result_files/varlo.parquet",
        "results/ceta_benchmark/Illumina_pe/called/ceta_multi/result_files/varlo.parquet",
        "results/ceta_benchmark/Illumina_pe/called/UNTREATED/result_files/varlo.parquet",
    output:
        "results/ceta_benchmark/Illumina_pe/called/UNTREATED/result_files/combined.txt",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/plots/plot_ceta_probs.log",
    shell:
        "touch {output}"