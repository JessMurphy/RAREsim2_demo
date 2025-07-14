#!/bin/bash

#SBATCH --output=%x_%a.out.%A
#SBATCH --error=%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=100-4000:100

# exit pipeline if a step fails
set -e

#cd /data001/projects/murphjes/RAREsim2_demo

# Define the population array
pop_list=(AFR EAS NFE SAS)

# Calculate the population index
pop_index=$(( ($SLURM_ARRAY_TASK_ID - 1) / 100 ))  # 0-3 for populations

# Access and export the variables using the task ID
export pop=${pop_list[$pop_index]}
export end=$(( $SLURM_ARRAY_TASK_ID - (100 * $pop_index) ))

echo "Task ID = $SLURM_ARRAY_TASK_ID, Pop: $pop,  Batch: $end"

TIME0=$(date +%s)

# a. Run RAREsim2 to create datasets for the same direction of effect and opposite directions of effect (50%/50%) scenarios
singularity exec "$container" ./3a_run_RAREsim2_power.sh

TIME1=$(date +%s)
echo "It took $(($TIME1 - $TIME0)) seconds to create datasets using RAREsim2 for the same direction of effect and opposite direction of effect scenarios (batches of 100 simulation replicates per population)."

# b. Run RAREsim2 to create datasets for the opposite directions of effect scenario with unequal amounts of risk/protective variants
singularity exec "$container" ./3b_run_RAREsim2_power_unequal.sh

TIME2=$(date +%s)
echo "It took $(($TIME2 - $TIME1)) additional seconds to create datasets using RAREsim2 for the opposite direction of effect scenario with unequal amounts of risk/protective variants (batches of 100 simulation replicates per population)."

# c. Run the rare variant association methods for all of the scenarios
singularity shell "$container" << EOF

TIME3=$(date +%s)

Rscript ./3c_run_methods_same_power.R "$SLURM_ARRAY_TASK_ID"

TIME4=$(date +%s)
echo "It took $(($TIME4 - $TIME3)) seconds to run the association methods for the same direction of effect scenario (batches of 100 simulation replicates per population)."

Rscript ./3c_run_methods_opp_power.R "$SLURM_ARRAY_TASK_ID"

TIME5=$(date +%s)
echo "It took $(($TIME5 - $TIME4)) seconds to run the association methods for the opposite directions of effect scenario (batches of 100 simulation replicates per population)."

Rscript ./3c_run_methods_opp_power_unequal.R "$SLURM_ARRAY_TASK_ID"

TIME6=$(date +%s)
echo "It took $(($TIME6 - $TIME5)) seconds to run the association methods for the opposite directions of effect scenario with unequal amounts of risk/protective variants (batches of 100 simulation replicates per population)."

EOF

TIME7=$(date +%s)
echo "It took $(($TIME7 - $TIME2)) seconds to run the association methods across the scenarios (batches of 100 simulation replicates per population)."


