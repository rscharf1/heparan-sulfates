#!/bin/bash
#SBATCH --job-name=homer_motif
#SBATCH --output=logs/homer_motif_%j.out
#SBATCH --error=logs/homer_motif_%j.err
#SBATCH --time=3:00:00              # Adjust if needed
#SBATCH --mem=64G                    # Memory allocation
#SBATCH --cpus-per-task=16           # Use 16 CPUs
#SBATCH --nodes=1
#SBATCH --ntasks=1

# findMotifsGenome.pl outputs/HS_promoters_mouse.bed mm10 outputs/homer -size given -mask

findMotifsGenome.pl inputs/cholesterol_promoters_human.bed hg38 outputs/homer_cholesterol_human -size given -mask -p 16 -preparsedDir outputs/homer_preparsed

findMotifsGenome.pl inputs/cholesterol_promoters_human.bed hg38 outputs/homer_cholesterol_human \
-size given -mask -genome ~/Tools/references/homer_genomes/hg38
