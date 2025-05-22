#!/bin/bash

#SBATCH --job-name=russ_job         # Job name
#SBATCH --output=slurm_log_%j.txt      # Standard output and error log
#SBATCH --mem=16G
#SBATCH --time=01:00:00             # Time limit hrs:min:sec

# Load R module (adjust to your environment)
module load R/4.4.2

# Run the R script
Rscript scripts/analysis.R
