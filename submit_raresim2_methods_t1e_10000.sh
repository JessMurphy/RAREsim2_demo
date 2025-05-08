#!/bin/bash

#SBATCH --output=/data001/projects/murphjes/RAREsim2/%x_%a.out.%A
#SBATCH --error=/data001/projects/murphjes/RAREsim2/%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=1000-40000:1000

# Debugging output
echo "Task ID: \"$SLURM_ARRAY_TASK_ID\""

singularity shell /storage/singularity/mixtures.sif << EOF

Rscript /data001/projects/murphjes/RAREsim2/raresim2_methods_t1e_10000.R "$SLURM_ARRAY_TASK_ID"


EOF