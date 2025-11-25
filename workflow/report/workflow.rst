Methylation Calling Evaluation with Varlociraptor
=================================================

This report summarizes the results of the *varlociraptor-methylation-evaluation* workflow, a Snakemake-based pipeline for benchmarking methylation calling on HG002 chromosome 21 (GRCh38, Ensembl 110).  
The workflow integrates methylation data from Illumina (bisulfite/EM-seq), Oxford Nanopore, and PacBio HiFi sequencing and compares Varlociraptor’s methylation calls against established platform-specific tools, namely Bismark, BSMAPz, MethylDackel, BisSNP, modkit, and pb-CpG-tools.

The generated metrics and visualizations provide a comprehensive assessment of methylation calling accuracy across technologies and callers, enabling reproducible, cross-platform evaluation of Varlociraptor’s performance.