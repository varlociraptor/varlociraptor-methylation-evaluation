# TODO: Do I have to hardcode the inputs, since the shell command uses each of them individually? 
rule call_methylation_together:
    input:
        varlo="resources/tools/varlociraptor/target/debug/varlociraptor",
        emseq1="results/Illumina_pe/EMSeq_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        emseq2="results/Illumina_pe/EMSeq_HG002_LAB02_{replicate}/normal_{scatteritem}.bcf",
        methylseq="results/Illumina_pe/MethylSeq_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        splat="results/Illumina_pe/SPLAT_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        truemethylbs="results/Illumina_pe/TrueMethylBS_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        truemethylox="results/Illumina_pe/TrueMethylOX_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        trueseq="results/Illumina_pe/TruSeq_HG002_LAB01_{replicate}/normal_{scatteritem}.bcf",
        scenario="resources/scenario_common.yaml",
    output:
        "results/common_calls/{replicate}/calls_{scatteritem}.bcf",
    log:
        "logs/varlociraptor/common_calls/{replicate}/call_methylation_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "{input.varlo} call variants generic --scenario {input.scenario} --obs emseq1={input.emseq1} emseq2={input.emseq2} methylseq={input.methylseq} splat={input.splat} truemethylbs={input.truemethylbs} truemethylox={input.truemethylox} trueseq={input.trueseq} > {output} 2> {log}"
