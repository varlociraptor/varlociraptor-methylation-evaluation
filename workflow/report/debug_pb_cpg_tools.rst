Left:
Strategy for identifying loci with high or low methylation that are not classified as fully methylated or unmethylated. We only use loci at which pb-CpG-tools reported a very low (< 5, top) or very high (> 95, bottom) methylation rate in at least one replicate. The shape indicates the classification: loci are labeled 'extreme' if Varlociraptor's methylation value equals $0$ (top) or equals $100$ (bottom), and 'unextreme' if it is at least $10$ (top) or is lower than $95$ (bottom). Loci with intermediate Varlociraptor values are excluded, as they do not clearly support either classification. Red dashed lines mark the thresholds used to select candidate loci. For all callers and platforms, we only use replicate values with a coverage of at least 20, computed with mosdepth independently of the callers.

Right:
Relative frequency
of methylation rates at the loci selected in A across all Illumina and Nanopore samples. To be independent of Varlociraptor, these rates are called by MethylDackel (Illumina) and Modkit (Nanopore).
The x-axis shows methylation level binned in increments of 5pp, the y-axis shows the relative proportion of calls in each bin. Colors correspond to the ‘low’/‘high’ classification derived from the Varlociraptor call on the PacBio sample.
