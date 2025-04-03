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

#########################################################
# Mason
#########################################################


# rule fake_meth_fasta_mason:
#     input:
#         "resources/chromosome_{chrom}.fasta",
#     output:
#         methylation="resources/Illumina_pe/simulated_data/chromosome_{chrom}_meth.fa",
#         variants="resources/Illumina_pe/simulated_data/chromosome_{chrom}_var.fa",
#         vcf="resources/Illumina_pe/simulated_data/chromosome_{chrom}.vcf",
#     conda:
#         "../envs/mason.yaml"
#     params:
#         pipeline_path=config["pipeline_path"],
#         # meth=lambda wildcards: 0.8 if wildcards.platform == "meth" else 1.0,
#         # snp=lambda wildcards: 0.0 if wildcards.platform == "meth" else 0.8,
#     shell:
#         """
#         mason_variator --in-reference {input} \
#                --methylation-levels \
#                --meth-cg-sigma 0.3 \
#                --meth-cg-mu 0.5 \
#                --meth-fasta-out {output.methylation} \
#                --out-fasta {output.variants} \
#                --out-vcf {output.vcf}

#         """


rule download_mason_latest:
    output:
        mason_dir=directory("resources/tools/seqan/apps/mason2"),
        mason="resources/tools/seqan/apps/mason2/methylation_levels.h",
    log:
        "../logs/download_mason.log",
    params:
        pipeline_path=config["pipeline_path"],
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p resources/tools
        cd resources/tools
        git clone git@github.com:seqan/seqan.git
        """


rule fake_methylation_mason:
    input:
        # mason=directory("resources/tools/seqan/apps/mason2"),
        chrom="resources/chromosome_{chrom}.fasta",
    output:
        methylation="resources/Illumina_pe/simulated_data/chromosome_{chrom}_meth.fa",
    conda:
        "../envs/mason.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        # meth=lambda wildcards: 0.8 if wildcards.platform == "meth" else 1.0,
        # snp=lambda wildcards: 0.0 if wildcards.platform == "meth" else 0.8,
    shell:
        """
        mason_methylation --in {input.chrom} \
            --methylation-levels \
            --meth-cg-sigma 0.3 \
            --meth-cg-mu 0.5 \
            --out {output.methylation} \
        """
        # cd {input.mason} && \


rule fake_variants_mason:
    input:
        # mason=directory("resources/tools/seqan/apps/mason2"),
        chrom="resources/chromosome_{chrom}.fasta",
    output:
        "resources/Illumina_pe/simulated_data/chromosome_{chrom}_variants.vcf",
    conda:
        "../envs/mason.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        # meth=lambda wildcards: 0.8 if wildcards.platform == "meth" else 1.0,
        # snp=lambda wildcards: 0.0 if wildcards.platform == "meth" else 0.8,
    shell:
        """
        mason_variator --in-reference {input.chrom} \
            --out-vcf {output} \
        """
        # cd {input.mason} && \
        # --snp-rate 0.01 \
        # --small-indel-rate 0.001 \
        # --sv-indel-rate 0.001 \


rule fake_reads_mason:
    input:
        genome="resources/chromosome_{chrom}.fasta",
        variants="resources/Illumina_pe/simulated_data/chromosome_{chrom}_variants.vcf",
        methylation="resources/Illumina_pe/simulated_data/chromosome_{chrom}_meth.fa",
    output:
        # bam="resources/Illumina_pe/simulated_data/alignment_{chrom}_test.bam",
        f1="resources/Illumina_pe/simulated_data/chromosome_{chrom}_f1.fastq",
        f2="resources/Illumina_pe/simulated_data/chromosome_{chrom}_f2.fastq",
        # truth="resources/Illumina_pe/simulated_data/chromosome_{chrom}.fasta",
    conda:
        "../envs/mason.yaml"
    params:
        pipeline_path=config["pipeline_path"],
        # meth=lambda wildcards: 0.8 if wildcards.platform == "meth" else 1.0,
        # snp=lambda wildcards: 0.0 if wildcards.platform == "meth" else 0.8,
    shell:
        """
        mason_simulator --input-reference {input.genome} \
                --num-fragments 10000 \
                --out {output.f1} \
                --out-right {output.f2} \
                --meth-fasta-in {input.methylation} \
                --enable-bs-seq \
                --illumina-read-length 150 \
        """
        # --seq-technology illumina \
        # --out-alignment {output.bam} \
        # --input-vcf {input.variants} \
        # mason_simulator --input-reference {input.variants} \


########Compute truth


# Mason has a different meth ratio for forward and reverse strands.
# That is why we need to compute the coverage on the forward and reverse strand independently.


rule align_simulated_reads_mason:
    input:
        fasta=expand(
            "resources/chromosome_{chrom}.fasta",
            chrom=config["simulated_chrom"],
        ),
        fasta_index=expand(
            "resources/chromosome_{chrom}.fasta.bwameth.c2t",
            chrom=config["simulated_chrom"],
        ),
        # fasta_index="resources/test_chr.fasta.bwameth.c2t",
        # fasta="resources/test_chr.fasta",
        f1=expand(
            "resources/Illumina_pe/simulated_data/chromosome_{chrom}_f1.fastq",
            chrom=config["simulated_chrom"],
        ),
        f2=expand(
            "resources/Illumina_pe/simulated_data/chromosome_{chrom}_f2.fastq",
            chrom=config["simulated_chrom"],
        ),
    output:
        "resources/Illumina_pe/simulated_data/alignment.sam",
    conda:
        "../envs/bwa-meth.yaml"
    threads: 30
    shell:
        """
        bwameth.py index-mem2 {input.fasta} && \
        bwameth.py --threads {threads} --reference {input.fasta} {input.f1} {input.f2} > {output}
        """


rule sam_bam_mason:
    input:
        "resources/Illumina_pe/simulated_data/alignment.sam",
    output:
        "resources/Illumina_pe/simulated_data/alignment.bam",
    conda:
        "../envs/samtools.yaml"
    threads: 30
    shell:
        """
        samtools view -Sb {input} > {output}
        """


rule mason_alignment_forward:
    input:
        "resources/Illumina_pe/simulated_data/alignment.bam",
    output:
        first="resources/Illumina_pe/simulated_data/alignment_99.bam",
        second="resources/Illumina_pe/simulated_data/alignment_147.bam",
        forward="resources/Illumina_pe/simulated_data/alignment_forward.bam",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -b -f 64 -F 16 {input} > {output.first}
        samtools view -b -f 16 -F 64 {input} > {output.second}
        samtools merge {output.forward} {output.first} {output.second}
        """


rule mason_alignment_reverse:
    input:
        "resources/Illumina_pe/simulated_data/alignment.bam",
    output:
        first="resources/Illumina_pe/simulated_data/alignment_83.bam",
        second="resources/Illumina_pe/simulated_data/alignment_163.bam",
        rev="resources/Illumina_pe/simulated_data/alignment_reverse.bam",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -b -f 16 -f 64 {input} > {output.first}
        samtools view -b -F 16 -F 64 {input} > {output.second}
        samtools merge {output.rev} {output.first} {output.second}
        """


rule sort_mason_reads:
    input:
        "resources/Illumina_pe/simulated_data/alignment_{orientation}.bam",
    output:
        # Name it like that in order to skip filtering on qual, mark_duplicates, ...
        "resources/Illumina_pe/simulated_data/alignment_sorted_{orientation}.bam",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    shell:
        """
        samtools sort -@ {threads}  {input} -o {output}    
        """


rule index_mason_alignment:
    input:
        "resources/Illumina_pe/simulated_data/alignment_sorted_{orientation}.bam",
    output:
        "resources/Illumina_pe/simulated_data/alignment_sorted_{orientation}.bam.bai",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"


rule mosdepth_mason_truth:
    input:
        bam="resources/Illumina_pe/simulated_data/alignment_sorted_{orientation}.bam",
        bai="resources/Illumina_pe/simulated_data/alignment_sorted_{orientation}.bam.bai",
        bed=lambda wildcards: expand(
            "resources/{chrom}/candidates.bed",
            chrom=config["simulated_chrom"],
        ),
    output:
        "resources/Illumina_pe/simulated_data/{orientation}_cov.mosdepth.global.dist.txt",
        "resources/Illumina_pe/simulated_data/{orientation}_cov.mosdepth.region.dist.txt",
        "resources/Illumina_pe/simulated_data/{orientation}_cov.regions.bed.gz",
        summary="resources/Illumina_pe/simulated_data/{orientation}_cov.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth_bed/{orientation}.log",
    params:
        extra="--no-per-base --use-median",  # optional
    # additional decompression threads through `--threads`
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v5.5.2/bio/mosdepth"


rule unzip_mosdepth_mason:
    input:
        "resources/Illumina_pe/simulated_data/{orientation}_cov.regions.bed.gz",
    output:
        "resources/Illumina_pe/simulated_data/{orientation}_cov.regions.bed",
    shell:
        """
        gunzip {input}
        """


rule mason_truth:
    input:
        cov_forward="resources/Illumina_pe/simulated_data/forward_cov.regions.bed",  # this named output is required for prefix parsing
        cov_reverse="resources/Illumina_pe/simulated_data/reverse_cov.regions.bed",  # this named output is required for prefix parsing
        methylation="resources/Illumina_pe/simulated_data/chromosome_{chrom}_meth.fa",
        candidates="resources/{chrom}/candidates.vcf",
    output:
        "resources/Illumina_pe/simulated_data/chromosome_{chrom}_truth.bed",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/ascii_to_meth.py"


###############################################################


# mason_simulator  \
#                 --input-reference {input.genome} \
#                 --meth-fasta-in {input.meth} \
#                 --num-fragments 10000 \
#                 --out {output.f1} \
#                 --out-right {output.f2} \
#                 --out-alignment {output.bam} \
#                 --seq-technology illumina \
#                 --methylation-levels \
#                 --enable-bs-seq \


###########################################################3
# --meth-fasta-out {output.truth} \
###########################################################3


#########################################################
# MethyFASTQ
#########################################################


# rule download_methylFastq:
#     output:
#         directory("resources/tools/MethylFASTQ"),
#     log:
#         "../logs/download_methylFastq.log",
#     params:
#         pipeline_path=config["pipeline_path"],
#     conda:
#         "../envs/install_program.yaml"
#     shell:
#         """
#         mkdir -p resources/tools
#         cd /resources/tools
#         git clone git@github.com:qBioTurin/MethylFASTQ.git
#         """


# # chg and chh = 1 because we are not interested in these anyways and want them to be unchanged in the bisulfite bams
# rule fake_meth_data:
#     input:
#         methylFastQ="resources/tools/MethylFASTQ/src/methylFASTQ.py",
#         # fasta="resources/test_chr.fasta",
#         fasta="resources/chromosome_{chrom}.fasta",
#     output:
#         f1="resources/Illumina_pe/simulated_data/chromosome_{chrom}_pe_f150r150_dir_R1.fastq",
#         f2="resources/Illumina_pe/simulated_data/chromosome_{chrom}_pe_f150r150_dir_R2.fastq",
#         ch3="resources/Illumina_pe/simulated_data/chromosome_{chrom}_pe_f150r150_dir.ch3",
#         # direc=directory("resources/Illumina_pe/simulated_data/"),
#     conda:
#         "../envs/methylFastQ.yaml"
#     params:
#         pipeline_path=config["pipeline_path"],
#         # meth=lambda wildcards: 0.8 if wildcards.platform == "meth" else 1.0,
#         # snp=lambda wildcards: 0.0 if wildcards.platform == "meth" else 0.8,
#     shell:
#         """
#         output_path=$(dirname "{output.f1}")
#         python {input.methylFastQ} -i {input.fasta} -o $output_path --seq paired_end --read 150 --chh 1.0 --chg 1.0 --cg 0.5 --snp 0.01 --error 0.005 --coverage 10
#         """
#         # python {input.methylFastQ} -i {input.fasta} -o /projects/koesterlab/benchmark-methylation/varlociraptor-methylation-evaluation/resources/Illumina_pe/simulated_data/{wildcards.sim_mode} --seq paired_end --read 4 --chh 1.0 --chg 1.0 --cg 1.0 --snp 0.0 --error 0.0 --coverage 2


###########################################################3
###########################################################3
# rule index_chromosome:
#     input:
#         "resources/chromosome_{chrom}.fasta",
#     output:
#         "resources/chromosome_{chrom}.fasta.bwameth.c2t",
#     conda:
#         "../envs/bwa-meth.yaml"
#     resources:
#         mem_mb=512,
#         runtime=10,
#     shell:
#         "bwameth.py index {input}"


# rule align_simulated_reads:
#     input:
#         fasta=expand(
#             "resources/chromosome_{chrom}.fasta",
#             chrom=config["simulated_chrom"],
#         ),
#         fasta_index=expand(
#             "resources/chromosome_{chrom}.fasta.bwameth.c2t",
#             chrom=config["simulated_chrom"],
#         ),
#         # fasta_index="resources/test_chr.fasta.bwameth.c2t",
#         # fasta="resources/test_chr.fasta",
#         f1=expand(
#             "resources/Illumina_pe/simulated_data/chromosome_{chrom}_pe_f150r150_dir_R1.fastq",
#             chrom=config["simulated_chrom"],
#         ),
#         f2=expand(
#             "resources/Illumina_pe/simulated_data/chromosome_{chrom}_pe_f150r150_dir_R2.fastq",
#             chrom=config["simulated_chrom"],
#         ),
#     output:
#         "resources/Illumina_pe/simulated_data/no_sra/alignment.bam",
#     conda:
#         "../envs/bwa-meth.yaml"
#     threads: 30
#     shell:
#         """
#         bwameth.py index-mem2 {input.fasta} && \
#         bwameth.py --threads {threads} --reference {input.fasta} {input.f1} {input.f2}  | samtools view -S -b - > {output}
#         """


rule sort_simulated_reads:
    input:
        "resources/Illumina_pe/simulated_data/alignment.bam",
    output:
        # Name it like that in order to skip filtering on qual, mark_duplicates, ...
        "resources/Illumina_pe/simulated_data/alignment_focused_downsampled_dedup_renamed.bam",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    shell:
        """
        samtools sort -@ {threads}  {input} -o {output}    
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
#         cargo run --release -- preprocess variants {input.chromosome} --candidates {input.candidates} --bam {input.alignments} --read-type $PLATFORM > {output}
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
#         cargo run --release -- call variants --omit-strand-bias generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output}
#         """
