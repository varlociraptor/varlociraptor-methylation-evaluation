# ref_tools = config["ref_tools"]


# # rule collect_ref_tools_illumina:
# #     input:
# #         bssnp="results/ref_tools/BisSNP/{protocol}/cpg.raw.CG.bedgraph",
# #         bsmap="results/ref_tools/bsMap/{protocol}/methylation_ratios.bed",
# #         bismark="results/ref_tools/bismark/{protocol}/alignment_bismark_sorted.bedGraph",
# #         dackel="results/ref_tools/dackel/{protocol}/alignments_CpG.bedGraph",
# #     output:
# #         tools_dir=directory("results/Illumina_pe/{protocol}/ref_methods"),
# #         bssnp="results/Illumina_pe/{protocol}/ref_methods/cpg.raw.CG.bedgraph",
# #         bsmap="results/Illumina_pe/{protocol}/ref_methods/methylation_ratios.bed",
# #         bismark="results/Illumina_pe/{protocol}/ref_methods/alignment_bismark_sorted.bedGraph",
# #         dackel="results/Illumina_pe/{protocol}/ref_methods/alignments_CpG.bedGraph",
# #     shell:
# #         """
# #         mkdir -p {output.tools_dir}
# #         cp {input.bssnp} {output.bssnp}
# #         cp {input.bsmap} {output.bsmap}
# #         cp {input.bismark} {output.bismark}
# #         cp {input.dackel} {output.dackel}
# #         """


# # rule ref_tools_illumina_finished:
# #     input:
# #         "results/ref_tools/BisSNP/{protocol}/cpg.raw.CG.bedgraph",
# #         "results/ref_tools/bsMap/{protocol}/methylation_ratios.bed",
# #         "results/ref_tools/bismark/{protocol}/alignment_bismark_sorted.bedGraph",
# #         "results/ref_tools/dackel/{protocol}/alignments_CpG.bedGraph",
# #     output:
# #         "results/Illumina_pe/{protocol}/ref_tools_completed.txt",
# #     shell:
# #         "touch {output}"


# rule collect_ref_tools_illumina:
#     input:
#         expand(
#             "results/ref_tools/{method}/{protocol}/{method}.bed",
#             protocol=config["data"]["Illumina_pe"].keys(),
#             method=config["ref_tools"]["Illumina_pe"],
#         ),
#     output:
#         directory(
#             expand(
#                 "results/Illumina_pe/{protocol}/ref_methods",
#                 protocol=config["data"]["Illumina_pe"].keys(),
#             )
#         ),
#         expand(
#             "results/Illumina_pe/{protocol}/ref_methods/{method}.bed",
#             protocol=config["data"]["Illumina_pe"].keys(),
#             method=config["ref_tools"]["Illumina_pe"],
#         ),
#     shell:
#         """
#         mkdir -p {output[0]}
#         cp {input} {output[1]}
#         """


# rule ref_tools_illumina_finished:
#     input:
#         expand(
#             "results/Illumina_pe/{protocol}/ref_methods/{method}.bed",
#             protocol=config["data"]["Illumina_pe"].keys(),
#             method=ref_tools["Illumina_pe"],
#         ),
#     output:
#         "results/Illumina_pe/{protocol}/ref_tools_completed.txt",
#     shell:
#         "touch {output}"


# rule collect_ref_tools_pacbio:
#     input:
#         expand(
#             "results/ref_tools/{method}/{protocol}/{method}.bed",
#             protocol=config["data"]["PacBio"].keys(),
#             method=config["ref_tools"]["PacBio"],
#         ),
#     output:
#         directory(
#             expand(
#                 "results/PacBio/{protocol}/ref_methods",
#                 protocol=config["data"]["PacBio"].keys(),
#             )
#         ),
#         expand(
#             "results/PacBio/{protocol}/ref_methods/{method}.bed",
#             protocol=config["data"]["PacBio"].keys(),
#             method=config["ref_tools"]["PacBio"],
#         ),
#     shell:
#         """
#         mkdir -p {output[0]}
#         cp {input} {output[1]}
#         """


# rule ref_tools_pacbio_finished:
#     input:
#         expand(
#             "results/PacBio/{protocol}/ref_methods/{method}.bed",
#             protocol=config["data"]["PacBio"].keys(),
#             method=ref_tools["PacBio"],
#         ),
#     output:
#         "results/PacBio/{protocol}/ref_tools_completed.txt",
#     shell:
#         "touch {output}"


# rule collect_ref_nanopore_illumina:
#     input:
#         expand(
#             "results/ref_tools/{method}/{protocol}/{method}.bed",
#             protocol=config["data"]["Nanopore"].keys(),
#             method=config["ref_tools"]["Nanopore"],
#         ),
#     output:
#         directory(
#             expand(
#                 "results/Nanopore/{protocol}/ref_methods",
#                 protocol=config["data"]["Nanopore"].keys(),
#             )
#         ),
#         expand(
#             "results/Nanopore/{protocol}/ref_methods/{method}.bed",
#             protocol=config["data"]["Nanopore"].keys(),
#             method=config["ref_tools"]["Nanopore"],
#         ),
#     shell:
#         """
#         mkdir -p {output[0]}
#         cp {input} {output[1]}
#         """
# rule ref_tools_nanopore_finished:
#     input:
#         expand(
#             "results/Nanopore/{protocol}/ref_methods/{method}.bed",
#             protocol=config["data"]["Nanopore"].keys(),
#             method=ref_tools["Nanopore"],
#         ),
#     output:
#         "results/Nanopore/{protocol}/ref_tools_completed.txt",
#     shell:
#         "touch {output}"
