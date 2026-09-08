#!/bin/bash

#SBATCH --job-name=bdg2bw
#SBATCH --output=bdg2bw_%j.out
#SBATCH --error=bdg2bw_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# Exit on error
set -e

# Loop through all .bedGraph files in the current directory
for file in inputs/*.bdg.txt; do
    echo "Processing $file"

    # Get the base filename (no extension)
    base=$(basename "$file" .bedGraph)

    echo "$base"

    sort -k1,1 -k2,2n inputs/$base > inputs/$base.sorted

    bedGraphToBigWig inputs/$base.sorted outputs/hg38.chrom.sizes outputs/$base.bw

    rm inputs/$base.sorted

    echo "Finished $file → outputs/${base}.bw"
done

echo "All files processed."