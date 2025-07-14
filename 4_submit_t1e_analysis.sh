#!/bin/bash

#SBATCH --output=%x_%a.out.%A
#SBATCH --error=%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=100-40000:100

# exit pipeline if a step fails
set -e

#cd /data001/projects/murphjes/RAREsim2_demo

# Define the population array
pop_list=(AFR EAS NFE SAS)

# Calculate the population index
pop_index=$(( ($SLURM_ARRAY_TASK_ID - 1) / 10000 ))  # 0-3 for populations

# Access and export the variables using the task ID
export pop=${pop_list[$pop_index]}
export end=$(( $SLURM_ARRAY_TASK_ID - (10000 * $pop_index) ))

echo "Task ID = $SLURM_ARRAY_TASK_ID, Pop: $pop,  Batch: $end"

# a. Run RAREsim2 to create additional datasets for the type I error scenario
if [ "$end" != "100" ]; then

TIME0=$(date +%s)

singularity exec "$container" ./4a_run_RAREsim2_t1e.sh

TIME1=$(date +%s)
echo "It took $(($TIME1 - $TIME0)) seconds to create additional type I error datasets with RAREsim2 (batches of 100 simulation replicates per population)."

fi

TIME2=$(date +%s)

# b. Run the rare variant association methods for the type I error scenario
singularity shell "$container" << EOF

Rscript ./4b_run_methods_t1e.R "$SLURM_ARRAY_TASK_ID"

EOF

TIME3=$(date +%s)
echo "It took $(($TIME3 - $TIME2)) seconds to run the association methods for the type I error scenario (batches of 100 simulation replicates per population)."

