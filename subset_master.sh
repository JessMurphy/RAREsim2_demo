#!/bin/bash

#SBATCH --job-name=subset_master_10000
#SBATCH --output=/data001/projects/murphjes/code/subset_master_10000.out.%j
#SBATCH --error=/data001/projects/murphjes/code/subset_master_10000.err.%j
#SBATCH -p math-alderaan
#SBATCH -N 1

singularity shell /storage/singularity/mixtures.sif << EOF 

Rscript /data001/projects/murphjes/code/subset_master.R

EOF

