rule download_mason2:
    output:
        "resources/tools/mason2-2.0.9-Linux-x86_64/bin/mason_methylation",
    log:
        "../logs/download_mason2.log",
    params:
        url="http://packages.seqan.de/mason2/mason2-2.0.9-Linux-x86_64.tar.xz",
        dest_dir="resources/tools",
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        mkdir -p {params.dest_dir}
        wget -O {params.dest_dir}/mason2-2.0.9-Linux-x86_64.tar.xz {params.url}
        tar -xf {params.dest_dir}/mason2-2.0.9-Linux-x86_64.tar.xz -C {params.dest_dir}
        """


rule run_mason_methylation_help:
    input:
        "resources/tools/mason2-2.0.9-Linux-x86_64/bin/mason_methylation",
    output:
        "test.txt",
    log:
        "../logs/mason_methylation_help_snake.log",
    conda:
        "../envs/install_program.yaml"
    shell:
        """
        {input} --help > {output}
        """


########################################################################################################################################3


rule fake_meth_data_methylFASTQ:
    input:
        fasta="resources/chromosome_1.fasta",
    output:
        "resources/example_paired/genome_pe_f150r150_dir_R1.fastq",
        "resources/example_paired/genome_pe_f150r150_dir_R2.fastq",
        "resources/example_paired/genome_pe_f150r150_dir.ch3",
    log:
        "logs/fake_meth_data.log",
    params:
        pipeline_path=config["pipeline_path"],
        methylFastQ_path=config["methylFastQ_path"],
    conda:
        "envs/methylFastQ.yaml"
    shell:
        """
        python {params.methylFastQ_path}src/methylFASTQ.py -i {input} -o {params.pipeline_path}resources/example_paired --seq paired_end --chh 1.0 --chg 1.0 --cg 0 --snp 0.0 --error 0.0 --coverage 4 
        """
