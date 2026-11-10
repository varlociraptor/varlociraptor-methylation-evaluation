rule build_varlociraptor:
    output:
        "resources/tools/meth_calling/varlociraptor/target/release/varlociraptor",
    log:
        "logs/build_varlociraptor.log",
    conda:
        "../envs/rust.yaml"
    params:
        branch=config.get("varlociraptor-branch"),
    shell:
        """
        (
        PARENT_DIR=$(dirname $(dirname $(dirname $(dirname {output}))))
        mkdir -p $PARENT_DIR
        cd $PARENT_DIR
        rm -rf varlociraptor
        git clone https://github.com/varlociraptor/varlociraptor.git
        cd varlociraptor
        git checkout {params.branch}
        export CC=x86_64-conda-linux-gnu-cc CXX=x86_64-conda-linux-gnu-c++
        cargo build --release
        ) &> {log}
        """


rule varlociraptor_preprocess:
    input:
        varlo=varlociraptor_binary(),
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
    params:
        varlo=varlociraptor_cmd,
    conda:
        "../envs/varlociraptor.yaml"
    resources:
        mem_mb=lambda wc: 6000 if wc.seq_platform in ["Nanopore", "PacBio"] else 2000,
    shell:
        """
        if [[ "{wildcards.seq_platform}" == "Illumina_pe" || "{wildcards.seq_platform}" == "Illumina_se" || "{wildcards.seq_platform}" == "Simulate" ]]; then
            meth="converted"
        else
            meth="annotated"
        fi
        {params.varlo} preprocess variants --omit-mapq-adjustment {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --methylation-read-type $meth --max-depth 5000 > {output} 2> {log}
        """


rule varlociraptor_call:
    input:
        varlo=varlociraptor_binary(),
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
        # Calls filtered by control_fdr_calls are not called again
        sample="(?!.*controlled_).*",
    params:
        varlo=varlociraptor_cmd,
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{params.varlo} call variants --no-absent generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output} 2> {log}"


# Filter the single- and multi-sample calls with varlociraptor's local FDR control on the bin events
# Local mode decides per record, so the scattered calls can be filtered separately before gathering.
rule control_fdr_calls:
    input:
        varlo=varlociraptor_binary(),
        calls="results/{call_type}/{seq_platform}/called/{sample}/calls_{scatteritem}.bcf",
    output:
        "results/{call_type}/{seq_platform}/called/{sample}/controlled_{fdr}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor_single/control_fdr_calls/{call_type}_{seq_platform}_{sample}_{fdr}_{scatteritem}.log",
    wildcard_constraints:
        sample="[^/]+",
        fdr=r"[0-9.]+",
    params:
        varlo=varlociraptor_cmd,
        events=" ".join(config["control_fdr"]["events"]),
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{params.varlo} filter-calls control-fdr {input.calls} --events {params.events} --fdr {wildcards.fdr} --smart-retain-artifacts --mode local-strict --subsets adjacent-5 > {output} 2> {log}"


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
