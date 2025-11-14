# Configuration Overview

## Genome Settings

Defines which genome build and version to use during alignment and reference file retrieval.

```yaml
sample:
  species: "homo_sapiens"
  datatype: "dna"
  build: "GRCh38"
  release: "110"
```

---

## Sequencing Platforms and Chromosome Selection

We focus our analysis on chromosome 21 for all sequencing platforms in order to avoid increasing the runtime and memory consumption. By commenting out a sequencing platform, the results for this platform will not be calculated.

```yaml
seq_platforms: 
  Illumina_pe: 21
  PacBio: 21
  Nanopore: 21
```

---

## Reference Methylation Callers for Benchmarking

* **Illumina** bisulfite data can be compared against: [bismark](https://doi.org/10.1093/bioinformatics/btr167), [bsMap](https://github.com/BSMAP/bsmap), [methylDackel](https://github.com/dpryan79/MethylDackel) and [bisSNP](https://doi.org/10.1186/gb-2012-13-7-r61)
* **PacBio** and **Nanopore** long-read platforms can be compared against: [modkit](https://github.com/PacificBiosciences/modkit), [pb_CpG_tools](https://github.com/PacificBiosciences/pb-CpG-tools)
* We evaluated Varlociraptor on multiple samples at the same time. To run the evaluation you have to define `multi_sample` and leave intentionally empty (required by downstream workflow)

```yaml
ref_tools:
  Illumina_pe: [bismark, bsMap, methylDackel]
  PacBio: [modkit, pb_CpG_tools]
  Nanopore: [modkit, pb_CpG_tools]
  multi_sample: []
```

---

## Input Data (Accessions and URLs)

We download the Illumina data using SRR accession numbers from the [EpiQC study](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02529-2). To keep the structure for every platform the same we use dummy names for PacBio and Nanopore.

```yaml
data:
  Illumina_pe:
    TruSeq_HG002_LAB01_REP01: [SRR13050995, SRR13050996, SRR13050998]
    TruSeq_HG002_LAB01_REP02: [SRR13050992, SRR13050993, SRR13050994]
    SPLAT_HG002_LAB01_REP01: [SRR13051055, SRR13051056, SRR13051057, SRR13051058]
    SPLAT_HG002_LAB01_REP02:  [SRR13051050, SRR13051051, SRR13051052, SRR13051054]
    MethylSeq_HG002_LAB01_REP01: [SRR13051104, SRR13051105, SRR13051106]
    MethylSeq_HG002_LAB01_REP02: [SRR13051101, SRR13051102, SRR13051103]
    EMSeq_HG002_LAB02_REP02: [SRR13051139]
    EMSeq_HG002_LAB02_REP01: [SRR13051140]
    EMSeq_HG002_LAB01_REP02: [SRR13051141]
    EMSeq_HG002_LAB01_REP01: [SRR13051142]
    TrueMethylOX_HG002_LAB01_REP02: [SRR13051230]
    TrueMethylOX_HG002_LAB01_REP01: [SRR13051231, SRR13051232]
    TrueMethylBS_HG002_LAB01_REP02: [SRR13051250, SRR13051251] 
    TrueMethylBS_HG002_LAB01_REP01: [SRR13051253] 
  PacBio:
    REP01: [pb_rep1]
    REP02: [pb_rep2]
  Nanopore:
    REP01: [np_rep1]
    REP02: [np_rep2]
  multi_sample:
    REP01: [dummy_rep1]
    REP02: [dummy_rep2]
```

The real PacBio and Nanopore data comes from direct download URLs:

```yaml
pb_rep1: https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam
pb_rep2: https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep2/analysis/HG002.m84005_220919_232112_s2.GRCh38.bam

np_rep1: https://42basepairs.com/download/s3/ont-open-data/giab_2025.01/basecalling/hac/HG002/PAW70337/calls.sorted.bam
np_rep2: https://42basepairs.com/download/s3/ont-open-data/giab_2025.01/basecalling/hac/HG002/PAW71238/calls.sorted.bam
```

---

## Sample Definitions

Names of the individual samples under which the workflow merges the two replicates.

```yaml
samples: 
  Illumina_pe: [TruSeq_HG002_LAB01, EMSeq_HG002_LAB02, EMSeq_HG002_LAB01, TrueMethylOX_HG002_LAB01, TrueMethylBS_HG002_LAB01, SPLAT_HG002_LAB01, MethylSeq_HG002_LAB01]
  PacBio: [REP]
  Nanopore: [REP]
  multi_sample: [REP]
```

---

## Methylation Calling Settings

* Filters reads below mapping quality 10.
* Enables parallelization using scatter-gather (20 tasks).
* Performs methylation calling at multiple FDR thresholds.

```yaml
min_mapping_quality: 10
scatter_number: 20
fdr_alpha: [0.01, 1.0]
```

---

## Plotting Settings

```yaml
heatmap_bin_size: 5
plot_type: svg
correlation_methods: [adjusted_mape]
```

* Binning resolution for heatmaps (bp).
* Output format for figures (`png`, `svg`, or `html`).
* Correlation metric(s) used to compare methylation calls.
