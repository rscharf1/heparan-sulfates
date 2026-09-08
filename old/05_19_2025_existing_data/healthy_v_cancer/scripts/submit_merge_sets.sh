#!/bin/bash

#SBATCH --job-name=russ_job         
#SBATCH --output=slurm_log_%j.txt
#SBATCH --mem=64G
#SBATCH --time=04:00:00

module load R/4.4.2

Rscript scripts/merge_sets.R