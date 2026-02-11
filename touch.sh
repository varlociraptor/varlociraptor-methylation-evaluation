BUCKET="s3://koesterlab/varlociraptor-methylation-evaluation/resources/Illumina_pe"


echo "Touching BAM files in $BUCKET"
rclone lsf -R "$BUCKET" | grep 'alignment_focused.bam$' | while read f; do
    echo "Touching $f"
    rclone touch "$BUCKET/$f"
done
