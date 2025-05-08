#!/bin/bash

#SBATCH --output=/data001/projects/murphjes/RAREsim2/%x_%a.out.%A
#SBATCH --error=/data001/projects/murphjes/RAREsim2/%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=0-3

# Define the arrays
pop_list=(AFR EAS NFE SAS)
NE_list=(17469 14269 11418 14269)
DL_list=(14627281 14673368 14705483 14508902) 

# Access and export the variables using the task ID
export pop=${pop_list[$SLURM_ARRAY_TASK_ID]}
export NE=${NE_list[$SLURM_ARRAY_TASK_ID]}
export DL=${DL_list[$SLURM_ARRAY_TASK_ID]}

echo "Task ID = $SLURM_ARRAY_TASK_ID, Pop: $pop"

singularity exec /storage/singularity/mixtures.sif /data001/projects/murphjes/RAREsim2/hapgen2_raresim2.sh