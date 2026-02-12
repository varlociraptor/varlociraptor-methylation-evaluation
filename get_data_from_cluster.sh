# Ich brauche genome.fasta, chromosome.fasta, 

rclone copy   ikim:/projects/koesterlab/benchmark-methylation/varlociraptor-methylation-evaluation/resources/Illumina_pe   s3:koesterlab/varlociraptor-methylation-evaluation/resources/Illumina_pe   --include "*/*/alignment_focused.bam"   --progress   -vv

rclone copy   ikim:/projects/koesterlab/benchmark-methylation/varlociraptor-methylation-evaluation/resources/genome.fasta   s3:koesterlab/varlociraptor-methylation-evaluation/resources/genome.fasta   --progress   -vv

rclone copy   ikim:/projects/koesterlab/benchmark-methylation/varlociraptor-methylation-evaluation/resources/chromosome_21.fasta   s3:koesterlab/varlociraptor-methylation-evaluation/resources/chromosome_21.fasta   --progress   -vv

rclone copy   ikim:/projects/koesterlab/benchmark-methylation/varlociraptor-methylation-evaluation/resources/ref_tools/bismark/bams   s3:koesterlab/varlociraptor-methylation-evaluation/resources/ref_tools/bismark/bams   --progress   -vv

rclone ls s3://koesterlab/varlociraptor-methylation-evaluation/resources/ref_tools/bismark/bams/ | \
awk '{print $2}' | \
while read f; do
    echo "Touching $f"
    rclone touch "s3://koesterlab/varlociraptor-methylation-evaluation/resources/ref_tools/bismark/bams/$f"
done

BUCKET="s3://koesterlab/varlociraptor-methylation-evaluation/resources/Illumina_pe"

# Liste alle alignment_focused.bam Dateien auf und touch sie
rclone lsf -R "$BUCKET" | grep 'alignment_focused.bam$' | while read f; do
    echo "Touching $f"
    rclone touch "$BUCKET/$f"
done
