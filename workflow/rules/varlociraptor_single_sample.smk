rule build_varlociraptor:
    output:
        "resources/tools/meth_calling/varlociraptor/target/release/varlociraptor",
    log:
        "logs/build_varlociraptor.log",
    conda:
        "../envs/rust.yaml"
    shell:
        """
        (
        PARENT_DIR=$(dirname $(dirname $(dirname $(dirname {output}))))
        mkdir -p $PARENT_DIR
        cd $PARENT_DIR
        rm -rf varlociraptor
        git clone https://github.com/varlociraptor/varlociraptor.git
        cd varlociraptor
        git checkout long-read-no-meth-info
        export CC=x86_64-conda-linux-gnu-cc CXX=x86_64-conda-linux-gnu-c++
        cargo build --release
        ) &> {log}
        """


rule varlociraptor_preprocess:
    input:
        chromosome=lambda wildcards: expand(
            "resources/{chrom}.fasta",
            chrom=chromosome_by_seq_platform.get(wildcards.seq_platform),
        ),
        genome_index=lambda wildcards: expand(
            "resources/{chrom}.fasta.fai",
            chrom=chromosome_by_seq_platform.get(wildcards.seq_platform),
        ),
        alignments="resources/{seq_platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam",
        alignment_index="resources/{seq_platform}/{sample}/candidate_specific/alignment_{scatteritem}.bam.bai",
        candidates=lambda wildcards: expand(
            "resources/{chrom}/candidates_{{scatteritem}}.bcf",
            chrom=(
                21
                if chromosome_by_seq_platform.get(wildcards.seq_platform) == "genome"
                else chromosome_by_seq_platform.get(wildcards.seq_platform)
            ),
        ),
    output:
        "results/preprocessed/{seq_platform}/{sample}/normal_{scatteritem}.bcf",
    log:
        "logs/varlociraptor_single/varlociraptor_preprocess/{seq_platform}_{sample}_{scatteritem}.log",
    benchmark:
        repeat(
            "benchmarks/{seq_platform}/varlociraptor/preprocessing/{sample}_{scatteritem}.bwa.benchmark.txt",
            config["benchmark_repeats"],
        )
    conda:
        "../envs/varlociraptor.yaml"
    resources:
        mem_mb=16000,
    shell:
        """
        if [[ "{wildcards.seq_platform}" == "Illumina_pe" || "{wildcards.seq_platform}" == "Illumina_se" || "{wildcards.seq_platform}" == "Simulate" ]]; then
            meth="converted"
        else
            meth="annotated"
        fi
        varlociraptor preprocess variants --omit-mapq-adjustment {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --methylation-read-type $meth --max-depth 5000 > {output} 2> {log}
        """


rule varlociraptor_call:
    input:
        preprocess_obs="results/preprocessed/{seq_platform}/{sample}/normal_{scatteritem}.bcf",
        scenario=workflow.source_path("../scenarios/scenario.yaml"),
    output:
        "results/single_sample/{seq_platform}/called/{sample}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor_single/varlociraptor_call/{seq_platform}_{sample}_{scatteritem}.log",
    benchmark:
        repeat(
            "benchmarks/{seq_platform}/varlociraptor/calling/{sample}_{scatteritem}.bwa.benchmark.txt",
            config["benchmark_repeats"],
        )
    wildcard_constraints:
        seq_platform="(?!multi_sample).*",
        sample="(?!simulated_data_no_bias).*",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor call variants generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output} 2> {log}"


rule gather_calls:
    input:
        gather.split_candidates(
            "results/{{call_type}}/{{seq_platform}}/{{sample}}/calls_{scatteritem}.vcf"
        ),
    output:
        "results/{call_type}/{seq_platform}/{sample}/result_files/varlo.bed",
    log:
        "logs/varlociraptor_single/gather_calls/{call_type}_{seq_platform}_{sample}.log",
    conda:
        "../envs/general.yaml"
    shell:
        "cat {input} > {output} 2> {log}"
