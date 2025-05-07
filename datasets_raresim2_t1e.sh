#!/bin/bash

nsim=10000 #number of individuals
num=19 #chromosome number
b=37 #block number
ncase=5000

WD=/data001/projects/murphjes/RAREsim2/datasets

start=$(($end-99))

# loop through the simulation replicates
for i in $(eval echo "{$start..$end}")
do 

echo "Simulation Replicate: ${i}"

# determine n from i
n=$(( (i - 1) / 1000 + 1 ))

prefix=${pop}/Round${n}/chr${num}.block${b}.${pop}.sim${i}

# convert the initial haplotype file to a sparse matrix
#python3 /home/math/murphjes/raresim2/convert.py \
#    -i ${WD}/Hapgen$((nsim/1000))K/${prefix}.controls.haps.gz \
#    -o ${WD}/Hapgen$((nsim/1000))K/${prefix}.controls.haps.sm

# prune functional and synonymous variants down to 100% functional / synonymous and remove the rows of zeros (use the z flag)
python3 /home/math/murphjes/raresim2/sim.py \
    -m ${WD}/Hapgen$((nsim/1000))K/${prefix}.controls.haps.sm \
    --functional_bins ${WD}/MAC_bins2/${pop}/Expected_MAC_${nsim}_${pop}_100_fun.txt \
    --synonymous_bins ${WD}/MAC_bins2/${pop}/Expected_MAC_${nsim}_${pop}_100_syn.txt \
    --stop_threshold 10 \
    -l ${WD}/Hapgen$((nsim/1000))K/${prefix}.copy.legend \
    -L ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.100fun.100syn.legend \
    -H ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.100fun.100syn.haps.gz \
    -z

# extract the type I error cases
python3 /home/math/murphjes/raresim2/extract.py \
    -i ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.100fun.100syn.haps.gz \
    -o ${WD}/Cases2/${prefix}.${ncase}.t1e.cases.100fun.100syn.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

# move the controls
mv ${WD}/Cases2/${prefix}.${ncase}.t1e.cases.100fun.100syn.haps-remainder.gz ${WD}/Controls2/${prefix}.${ncase}.controls.same.100fun.100syn.haps-remainder.gz

done

# type I error cases: chr19.block37.${pop}.sim${rep}.${ncase}.t1e.cases.100fun.100syn.haps-sample.gz
# controls (same): chr19.block37.${pop}.sim${rep}.${ncase}.t1e.cases.100fun.100syn.haps-remainder.gz