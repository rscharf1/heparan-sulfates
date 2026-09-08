#!/bin/bash

#SBATCH --job-name=PCHi-C
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --time=3:00:00              # Adjust if needed
#SBATCH --mem=64G                    # Memory allocation
#SBATCH --cpus-per-task=16           # Use 16 CPUs
#SBATCH --nodes=1
#SBATCH --ntasks=1

# bedtools intersect to find regions of enhancer/promoter interactions within active chromatin regions

mkdir tmp

# Filter interactions to have meaningful reads (>= 5)
awk 'BEGIN{OFS="\t"} $9 >= 5 {print $5, $6, $7, $4}' \
inputs/HPC7_Promoter_Capture_Interactions.ibed > tmp/filtered_otherEnds_with_gene.bed

tail -n +2 tmp/filtered_otherEnds_with_gene.bed > tmp/filtered_otherEnds_with_gene.bed.tmp 

mv tmp/filtered_otherEnds_with_gene.bed.tmp tmp/filtered_otherEnds_with_gene.bed

awk '{OFS="\t"; $1="chr"$1; print}' tmp/filtered_otherEnds_with_gene.bed > tmp/filtered_otherEnds_with_gene.bed.tmp

mv tmp/filtered_otherEnds_with_gene.bed.tmp tmp/filtered_otherEnds_with_gene.bed

# Filter and merge H3K27ac 
awk 'BEGIN{OFS="\t"} $4 >= 100 {print $1, $2, $3, $4}' \
inputs/GSM1329815_ChIPseq_H3K27Ac_HPC7.cleaned.bedgraph > tmp/H3K27ac.bed

bedtools merge -i inputs/GSM1329815_ChIPseq_H3K27Ac_HPC7.cleaned.bedgraph -c 4,4 -o max,mean > tmp/H3K27ac_merged_with_signal.bed

# Merge datasets
bedtools intersect -a tmp/filtered_otherEnds_with_gene.bed \
                   -b tmp/H3K27ac_merged_with_signal.bed \
                   -wa -wb > outputs/active_enhancers_with_genes.bed
