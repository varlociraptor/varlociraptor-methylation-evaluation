Example of the formula tree that Varlociraptor traverses to compute the joint probability Pr(Z_s₁, Z_s₂ | e) Pr(e) of an event e over two samples s₁ and s₂, here for the event ``s1: ]0.0,1.0] & (s2:0.0 | l2fc(s1, s2) > 0.5)``. That means s₁ is methylated and s₂ is either unmethylated or at least 2^0.5 ≈ 1.41-fold less methylated than s₁.
Each node constrains the methylation rate of one sample, and the two branches correspond to the two alternatives of the disjunction (``|``).
The tree is traversed depth-first. The root assigns s₁ a methylation rate in ]0, 1], i.e. the recursion descends once for each grid point of the numerical integration over this interval.
On the left path, s₂ is assigned the discrete value 0.
On the right path, the log2 fold change does not assign a value itself but narrows the interval of s₂, given the value of s₁, to s₂ < s₁ / 2^0.5.
Since the formula does not constrain s₂ further on this path, s₂ is added according to its universe [0, 1] (bottom right), such that each path from root to leaf constrains every sample exactly once.
Only at the last node of a path, where each sample has a single methylation rate, the joint probability is evaluated as the product of the likelihoods and the prior defined by the scenario, and returned to the caller.
