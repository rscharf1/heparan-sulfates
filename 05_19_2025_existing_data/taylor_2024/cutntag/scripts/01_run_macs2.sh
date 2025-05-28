#!/bin/bash

# Exit on error
set -e

# Loop through all .bedGraph files in the current directory
for file in inputs/*.bedgraph; do
    echo "Processing $file"

    # Get the base filename (no extension)
    base=$(basename "$file" .bedGraph)

    echo "$base"

    # # Call peaks with MACS2
    macs2 bdgpeakcall \
        -i "inputs/${base}" \
        -c 2 -l 100 \
        -o "outputs/${base}_peaks.bed"
    
    echo "Finished $file → outputs/${base}_peaks.bed"
done

echo "All files processed."