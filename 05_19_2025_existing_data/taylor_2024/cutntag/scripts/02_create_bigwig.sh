#!/bin/bash

# Exit on error
set -e

# Loop through all .bedGraph files in the current directory
for file in outputs/*.bed; do
    echo "Processing $file"

    # Get the base filename (no extension)
    base=$(basename "$file" .bedGraph)

    echo "$base"

    # cut -f-4 outputs/${base} > outputs/${base}.signal

    awk 'BEGIN {OFS="\t"} !/^track/ {print $1, $2, $3, $5}' \
        outputs/${base} > outputs/${base}.signal

    sort -k1,1 -k2,2n outputs/${base}.signal > outputs/${base}.signal.sorted

    grep -v "^track" outputs/${base}.signal.sorted > outputs/${base}.signal.sorted.clean

    bedGraphToBigWig \
        "outputs/${base}.signal.sorted.clean" \
        inputs/hg38.chrom.sizes \
        "outputs/${base}_peaks.bw"

    rm outputs/${base}.signal
    rm outputs/${base}.signal.sorted
    rm outputs/${base}.signal.sorted.clean
    
    echo "Finished $file → outputs/${base}.bed"
done

echo "All files processed."