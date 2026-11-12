Calibration of Varlociraptor's posterior distribution of the methylation rate θ on simulated data.
For every CpG site, u = P(θ ≤ θ* | data) is computed from the posterior, where θ* is the true methylation rate.
The posterior consists of a point mass at θ = 0 (event ABSENT) and the continuous allele frequency distribution (AFD) on ]0,1].
The AFD is integrated with the trapezoidal rule between its grid points; its entry at θ = 0 is ignored, since that point mass is given by PROB_ABSENT.
For sites with θ* = 0, u is drawn uniformly from [0, P(θ = 0 | data)].
The histogram shows the number of sites per bin of u. The dashed line marks the count expected under perfect calibration (uniform u).
