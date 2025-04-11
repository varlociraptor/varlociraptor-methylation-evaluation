
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/varlociraptor/varlociraptor-methylation-evaluation/workflows/Tests/badge.svg?branch=main)](https://github.com/varlociraptor/varlociraptor-methylation-evaluation/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for the evaluation of the methylation calling of varlociraptor.

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=varlociraptor%2Fvarlociraptor-methylation-evaluation).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

# Methylation Calling with Varlociraptor

**varlociraptor-methylation-evaluation** is a Snakemake-based workflow designed for accurate methylation calling using [Varlociraptor](https://github.com/varlociraptor/varlociraptor). This pipeline integrates simulation, evaluation, and visualization tools for benchmarking methylation data from various sequencing technologies, with a strong emphasis on configuration flexibility and reproducibility.

---

## Overview

This workflow focuses on methylation rate estimation via Varlociraptor and compares it with other tools across multiple platforms such as:

- Illumina paired-end (PE)
- PacBio
- Nanopore

Key features include:

- Simulation-ready configuration
- Platform-specific reference tool comparison
- Coverage binning and statistical plots
- Configurable thresholds and methylation types
- Ready-to-use with GRCh38 (release 110)

---

## Configuration

Configuration is stored in a YAML format and defines the full experimental setup. Below is a breakdown of the key sections from your `config.yaml`.

### 🧬 Genome & Sample Info

```yaml
sample:
  species: "homo_sapiens"
  datatype: "dna"
  build: "GRCh38"
  release: "110"
