#!/bin/bash

#SBATCH --output=%x_%a.out.%A
#SBATCH --error=%x_%a.err.%A
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=100-40000:100

# exit pipeline if a step fails
set -e

# calculate necessary variables
pop_index=$(( ($SLURM_ARRAY_TASK_ID - 1) / 10000 ))  # 0-3 for populations
end=$(( $SLURM_ARRAY_TASK_ID - (10000 * $pop_index) ))

echo "Task ID = $SLURM_ARRAY_TASK_ID, Pop Index: $pop_index, Batch: $end"

# a. Make master legend file (only need to do once)
if [ "$end" == "100" ]; then

TIME0=$(date +%s)

singularity shell "$container" << EOF

Rscript 1a_make_master_legend.R "$SLURM_ARRAY_TASK_ID"

EOF

TIME1=$(date +%s)
echo "It took $(($TIME1 - $TIME0)) seconds to make a master legend file for one population."

fi

TIME2=$(date +%s)

# b. Subset master legend file
singularity shell "$container" << EOF

Rscript 1b_subset_master_legend.R "$SLURM_ARRAY_TASK_ID"

EOF

TIME3=$(date +%s)
echo "It took $(($TIME3 - $TIME2)) seconds to make 100 legend files for one population."