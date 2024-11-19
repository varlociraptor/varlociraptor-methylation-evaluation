# """
# Mein Plan:
# Erstelle 2 fake alignment files
#     1. File mit methylierung, snp = 0, methylierung Cpg=x
#     2. File mit snp, snp = x, methylierung= 1 (Damit alles C bleibt und von Varlo bei SNPs nicht als T aufgegriffen wird)
# Rufe die alignments mit Varlo auf
#     1. Mit normalen candidate File
#     2. Mit bearbeitetem Candidate file, anstelle von Meth muss da was auch immer stehen
# Erstelle Scenario File
# Kann ich Scatteritems miteinander reinpacken und nicht die Ganze? Sollte doch eigtl klappen, da alle gleich lang sind

# Wie interpretiere ich Ergebnis?
# """


rule download_methylFastq:
    output:
        directory("resources/tools/MethylFASTQ"),
    log:
        "../logs/download_methylFastq.log",
    params:
        pipeline_path=config["pipeline_path"],
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p {params.pipeline_path}resources/tools
        cd {params.pipeline_path}/resources/tools
        git clone git@github.com:qBioTurin/MethylFASTQ.git
        """


# chg and chh = 1 because we are not interested in these anyways and want them to be unchanged in the bisulfite bams
rule fake_meth_data:
    input:
        methylFastQ="resources/tools/MethylFASTQ/src/methylFASTQ.py",
        fasta="resources/test_chr.fasta",
    output:
        f1="resources/Illumina_pe/simulated_data/{sim_mode}/test_chr_pe_f4r4_dir_R1.fastq",
        f2="resources/Illumina_pe/simulated_data/{sim_mode}/test_chr_pe_f4r4_dir_R2.fastq",
        ch3="resources/Illumina_pe/simulated_data/{sim_mode}/test_chr_pe_f4r4_dir.ch3",
    conda:
        "../envs/methylFastQ.yaml"
    log:
        "logs/fake_meth_data{sim_mode}.log",
    params:
        pipeline_path=config["pipeline_path"],
        # meth=lambda wildcards: 0.8 if wildcards.platform == "meth" else 1.0,
        # snp=lambda wildcards: 0.0 if wildcards.platform == "meth" else 0.8,
    shell:
        """
        python {input.methylFastQ} -i {input.fasta} -o /projects/koesterlab/benchmark-methylation/varlociraptor-methylation-evaluation/resources/Illumina_pe/simulated_data/{wildcards.sim_mode} --seq paired_end --read 4 --chh 1.0 --chg 1.0 --cg 1.0 --snp 0.0 --error 0.0 --coverage 2
        """


rule index_genome:
    input:
        "resources/test_chr.fasta",
    output:
        "resources/test_chr.fasta.bwameth.c2t",
    log:
        "logs/index_genome.log",
    conda:
        "../envs/bwa-meth.yaml"
    shell:
        "bwameth.py index-mem2 {input}"


rule align_simulated_reads:
    input:
        fasta_index="resources/test_chr.fasta.bwameth.c2t",
        fasta="resources/test_chr.fasta",
        reads1="resources/Illumina_pe/simulated_data/{sim_mode}/test_chr_pe_f4r4_dir_R1.fastq",
        reads2="resources/Illumina_pe/simulated_data/{sim_mode}/test_chr_pe_f4r4_dir_R1.fastq",
    output:
        "resources/Illumina_pe/simulated_data/{sim_mode}/alignment.bam",
    conda:
        "../envs/bwa-meth.yaml"
    log:
        "logs/align_simulated_reads_{sim_mode}.log",
    threads: 30
    resources:
        mem_mb=512,
    shell:
        """
        bwameth.py --threads {threads} --reference {input.fasta} {input.reads1} {input.reads2}  | samtools view -S -b - > {output}
        """


# rule compute_meth_observations:
#     input:
#         chromosome="resources/chromosome_21.fasta",
#         genome_index="resources/chromosome_21.fasta.fai",
#         alignments="resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam",
#         alignment_index="resources/{platform}/{protocol}/candidate_specific/alignment_valid_{scatteritem}.bam.bai",
#         candidates=lambda wildcards: expand(
#             "resources/{chrom}/candidates_{{scatteritem}}.bcf",
#             chrom=chromosome_by_platform[wildcards.platform],
#         ),
#     output:
#         temp("results/{platform}/{protocol}/normal_{scatteritem}.bcf"),
#     log:
#         "logs/compute_meth_observations_{platform}_{protocol}_{scatteritem}.log",
#     conda:
#         "envs/varlociraptor.yaml"
#     params:
#         varlo_path=config["varlo_path"],
#         pipeline_path=config["pipeline_path"],
#     shell:
#         """
#         cd {params.varlo_path}
#         if [[ "{wildcards.platform}" == "Illumina_pe" || "{wildcards.platform}" == "Illumina_se" ]]; then
#             PLATFORM="Illumina"
#         else
#             PLATFORM="{wildcards.platform}"
#         fi
#         cargo run --release -- preprocess variants {params.pipeline_path}{input.chromosome} --candidates {params.pipeline_path}{input.candidates} --bam {params.pipeline_path}{input.alignments} --read-type $PLATFORM > {params.pipeline_path}{output}
#         """
# rule call_varlo_joint:
#     input:
#         preprocess_obs="results/Illumina_pe/fake_meth/normal_{scatteritem}.bcf",
#         preprocess_obs="results/Illumina_pe/fake_snp/normal_{scatteritem}.bcf",
#         scenario="resources/scenario_joint.yaml",
#     output:
#         temp("results/{platform}/{protocol}/calls_{scatteritem}.bcf"),
#     log:
#         "logs/call_methylation_{platform}_{protocol}_{scatteritem}.log",
#     conda:
#         "envs/varlociraptor.yaml"
#     params:
#         varlo_path=config["varlo_path"],
#         pipeline_path=config["pipeline_path"],
#     shell:
#         """
#         cd {params.varlo_path}
#         cargo run --release -- call variants --omit-strand-bias generic --scenario {params.pipeline_path}{input.scenario} --obs normal={params.pipeline_path}{input.preprocess_obs} > {params.pipeline_path}{output}
#         """
