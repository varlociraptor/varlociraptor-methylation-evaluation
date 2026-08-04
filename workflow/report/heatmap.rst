{% set platform_labels = {
    "np_pb": "Nanopore and PacBio",
    "np_methylSeq": "Nanopore and Illumina",
    "pb_methylSeq": "PacBio and Illumina",
    "Nanopore": "Nanopore",
    "PacBio": "PacBio",
    "Illumina_pe": "Illumina",
} %}



Discordance between predicted methylation rates of each pair of replicates {% if snakemake.wildcards.sample == "all_samples" %}across all Illumina protocols{% else %}on {{ platform_labels.get(snakemake.wildcards.seq_platform, snakemake.wildcards.seq_platform) }} data{% endif %}, stratified by caller. Each squared bin counts the number of prediction pairs within the associated methylation rate ranges on the two axes. The fewer counts at the top left and bottom right edges and the more on the diagonal, the better is the concordance of the callers predictions across the replicate pairs.

{% if snakemake.wildcards.call_type == "multi_sample" %}

We combined two single {{ platform_labels.get(snakemake.wildcards.seq_platform, snakemake.wildcards.seq_platform) }} samples as joint evidence for multi-sample methylation calling using Varlociraptor. By leveraging information across multiple samples, Varlociraptor increases confidence in the computed methylation rates, which is reflected in a sharper diagonal across a larger number of loci N.



{% endif %}
