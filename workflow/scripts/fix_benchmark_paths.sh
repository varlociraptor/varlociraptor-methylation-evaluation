#!/bin/bash

BENCHMARK_DIR="/home/adrian/Documents/varlociraptor-methylation-evaluation/benchmarks/Illumina_pe/bismark"
OUTPUT_DIR="$BENCHMARK_DIR/bismark_align"

mkdir -p "$OUTPUT_DIR"

# Step 1: Copy all files with SRR ID
for dir in "$BENCHMARK_DIR"/*_SRR*/; do
    [ -d "$dir" ] || continue
    srr_id=$(basename "$dir" | grep -o "SRR[0-9]*")
    for file in "$dir"/*.bwa.benchmark.txt; do
        [ -f "$file" ] || continue
        base=$(basename "$file" .bwa.benchmark.txt)
        cp "$file" "$OUTPUT_DIR/${base}_${srr_id}.bwa.benchmark.txt"
        echo "Copied: ${base}_${srr_id}.bwa.benchmark.txt"
    done
done

# Step 2: Get unique base names and concatenate
cd "$OUTPUT_DIR"
ls *_SRR*.bwa.benchmark.txt | sed 's/_SRR[0-9]*\.bwa\.benchmark\.txt$//' | sort -u > /tmp/bases.txt

while read base; do
    srr_files="${base}_SRR*.bwa.benchmark.txt"
    count=$(ls -1 $srr_files 2>/dev/null | wc -l)

    if [ "$count" -gt 1 ]; then
        echo "Concatenating $base.bwa.benchmark.txt ($count files)..."
        head -1 $(ls -1 $srr_files | head -1) > "$base.bwa.benchmark.txt"
        for f in $srr_files; do
            tail -n +2 "$f" >> "$base.bwa.benchmark.txt"
        done
        rm $srr_files
        echo "Done: $base.bwa.benchmark.txt"
    fi
done < /tmp/bases.txt

rm /tmp/bases.txt
echo "Finished!"
