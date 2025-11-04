
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A59.10.1-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/varlociraptor/varlociraptor-methylation-evaluation/workflows/Tests/badge.svg?branch=main)](https://github.com/varlociraptor/varlociraptor-methylation-evaluation/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for the evaluation of the methylation calling of varlociraptor.


# Methylation Calling with Varlociraptor

**varlociraptor-methylation-evaluation** is a Snakemake-based workflow designed for accurate methylation calling using [Varlociraptor](https://github.com/varlociraptor/varlociraptor).  The evaluation is performed on chromosome 21 of the human reference genome GRCh38 (Ensembl release 110) using methylation data generated from three sequencing technologies: Illumina (bisulfite/EM-seq sequencing), Oxford Nanopore, and PacBio HiFi.

To benchmark Varlociraptor, we compare its methylation calls against established, platform-specific methylation callers:

* Illumina sequencing

  * [Bismark](https://doi.org/10.1093/bioinformatics/btr167)
  * [BSMAP](https://github.com/BSMAP/bsmap)
  * [MethylDackel](https://github.com/dpryan79/MethylDackel)
  * [BisSNP](https://doi.org/10.1186/gb-2012-13-7-r61)

* PacBio and Oxford Nanopore sequencing

  * [modkit](https://github.com/PacificBiosciences/modkit)
  * [pb_CpG_tools](https://github.com/PacificBiosciences/pb-CpG-tools)

The workflow produces multiple benchmarking outputs, enabling a comprehensive comparison of methylation calling accuracy across sequencing platforms and tools.

---

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=varlociraptor%2Fvarlociraptor-methylation-evaluation).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

---

## Configuration

Configuration is stored in a YAML format and defines the full experimental setup. See this [README](config/README.md) for a detailed explanation of the configuration options.