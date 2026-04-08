#!/usr/bin/env bash
set -euo pipefail

SOURCE="ikim:/projects/koesterlab/benchmark-methylation/varlociraptor-methylation-evaluation-old"
DEST="s3:koesterlab/varlociraptor-methylation-evaluation"
PREFIX="resources/Illumina_pe"

FILES=(
    # "TrueMethylBS_HG002_LAB01_REP01/SRR13051253/SRR13051253"
    # "TrueMethylOX_HG002_LAB01_REP01/SRR13051232/SRR13051232"
    # "SPLAT_HG002_LAB01_REP01/SRR13051056/SRR13051056"
    # "SPLAT_HG002_LAB01_REP02/SRR13051050/SRR13051050"
    # "SPLAT_HG002_LAB01_REP02/SRR13051051/SRR13051051"
    "MethylSeq_HG002_LAB01_REP01/SRR13051104/SRR13051104"
    "EMSeq_HG002_LAB01_REP01/SRR13051142/SRR13051142"
)

for file in "${FILES[@]}"; do
    for i in 1 2; do
        echo $SOURCE/$PREFIX/${file}_${i}_trimmed.fastq \
            $DEST/$PREFIX/$(dirname "${file}")
        rclone copy \
            $SOURCE/$PREFIX/${file}_${i}_trimmed.fastq \
            $DEST/$PREFIX/$(dirname "${file}") \
            --progress --transfers 4
    done
done
