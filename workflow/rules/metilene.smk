rule download_metilene:
    output:
        temp("resources/tools/metilene_v02-8.tar.gz"),
    shell:
        """
        wget -O {output} http://www.bioinf.uni-leipzig.de/Software/metilene/metilene_v02-8.tar.gz
        """


rule unpack_metilene:
    input:
        "resources/tools/metilene_v02-8.tar.gz",
    output:
        directory("resources/tools/metilene"),
    shell:
        """
        mkdir -p {output}
        tar -xzf {input} -C {output} --strip-components=1
        """


rule metilene_input:
    input:
        "results/Illumina_pe/TruSeq_HG002_LAB01_REP01/calls.vcf",
        "results/{group2}/data1/calls.vcf",
    output:
        "results/dmr_calls/metilene_input_{group2}.txt",
    log:
        "../logs/metilene_input_{group2}.log",
    script:
        "../scripts/metilene_input.py"


rule call_metilene:
    input:
        "results/dmr_calls/metilene_input_{group2}.txt",
    output:
        "results/dmr_calls/metilene_output_{group2}.txt",
    log:
        "../logs/call_metilene_{group2}.log",
    conda:
        "../envs/metilene.yaml"
    threads: 10
    shell:
        """
        metilene -t {threads} -a Illumina_pe -b {wildcards.group2} {input} > {output}
        """


rule plot_dmrs:
    input:
        met=directory("resources/tools/metilene"),
        met_out="results/dmr_calls/metilene_output_{group2}.txt",
    output:
        "results/dmr_calls/{group2}_qval.0.05.bedgraph",
        "results/dmr_calls/{group2}_qval.0.05.pdf",
    conda:
        "../envs/metilene.yaml"
    params:
        prefix=lambda wildcards: f"results/dmr_calls/{wildcards.group2}",
    shell:
        """
        perl {input.met}/metilene_output.pl -q {input.met_out} -o {params.prefix} -a Illumina_pe -b {wildcards.group2}
        """
