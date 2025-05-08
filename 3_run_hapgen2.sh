#!/bin/bash

# variables to be updated per simulation
Nsim=10000 #number of individuals
num=19 #chromosome number
b=37 #block number

# input variables
#WD=/data001/projects/murphjes
Hap=./input/1000G/${pop}_Block${b}_CDS_ref_added.hap
Map=./input/1000G/genetic_map_chr${num}_combined.txt

start=1
end=5

for i in $(eval echo "{$start..$end}")
do

# determine n from i
n=$(( (i - 1) / 1000 + 1 ))

# output variables
Output=./datasets/Hapgen$((Nsim/1000))K/${pop}/Round${n}/chr${num}.block${b}.${pop}.sim${i}

# save a copy of the legend file
cp $Output.legend $Output.copy.legend

# unzip the haplotype file
gunzip $Hap.gz

# simulate with HAPGEN2
/storage/singularity/mixtures.sif hapgen2 -h $Hap \
-m $Map \
-l $Output.copy.legend \
-o $Output.gz \
-n $Nsim 0 \
-dl $DL 1 1.0 1.0 \
-Ne $NE \
-no_gens_output

# remove unnecessary files produced by Hapgen
rm $Output.cases.*
rm $Output.controls.sample

echo $i

done