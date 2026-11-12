Number of CpG loci with a confident bias called by Varlociraptor, by sequencing technology.
For this plot, the calls of each replicate are filtered separately with Varlociraptor's FDR control on the ARTIFACT
event only (``varlociraptor filter-calls control-fdr --events ARTIFACT --mode local-strict``), at the level α given in
the row title. A replicate counts as biased at a locus if the locus passes this filter and
Varlociraptor reports a bias field (SB, ROB, RPB, SCB, HE, ALB or FEB) for it, i.e. the locus is confidently an
artifact at FDR α.
Panels show Illumina with all samples pooled (loci counted per sample), PacBio and Nanopore


Each locus is counted once. The bias type (color) is that of the confident replicate. If both are confident and the
bias types differ, the type of replicate 1 is shown. Only loci called in both replicates are considered.
