#!/bin/bash

#SBATCH --job-name=homer
#SBATCH --time=3:00:00              
#SBATCH --mem=64G                  
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

export HOMER_HOME=/gs/gsfs0/users/rscharf/Tools/homer
export PATH=$HOMER_HOME/bin:$PATH

~/Tools/homer/bin/findMotifs.pl inputs/lysosomal_genes.txt human outputs/lysosomal_genes/ -start -400 -end 100 -p8
