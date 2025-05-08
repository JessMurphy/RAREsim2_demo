#!/bin/bash

#singularity shell /storage/singularity/mixtures.sif

nsim=10000 #number of individuals
num=19 #chromosome number
b=37 #block number
ncase=5000 #number of cases (also number of controls)
pcase=(160 140 130 120 110) #percent of expected functional variants in cases
same=(160 140 120) #pcase necessary for same direction of effect scenario
opp=(130 120 110) #pcase necessary for opposite directions of effect scenario

WD=/data001/projects/murphjes/RAREsim2/datasets

# CREATE EXPECTED MAC BINS (need target data for all populations)

#for p in "${!pcase[@]}"; do

#/storage/singularity/mixtures.sif python3 /home/math/murphjes/raresim2/expected_vars.py \
#    --mac ${WD}/MAC_bins2/MAC_bins_${nsim}.txt \
#    -o ${WD}/MAC_bins2/${pop}/Expected_MAC_${nsim}_${pop}_${pcase[$p]}.txt \
#    -N $((2*$ncase)) \
#    --nvar_target_data ${WD}/MAC_bins2/${pop}/chr19_block37_${pop}_nvar_target_data.txt \
#    --afs_target_data ${WD}/MAC_bins2/${pop}/chr19_block37_${pop}_AFS_target_data.txt \
#    --w_fun $(echo "scale=1; ${pcase[$p]} / 100" | bc) \
#    --reg_size 19.029

#rm ${WD}/MAC_bins2/${pop}/Expected_MAC_${nsim}_${pop}_${pcase[$p]}_syn.txt

#done

#/storage/singularity/mixtures.sif python3 /home/math/murphjes/raresim2/expected_vars.py \
#    --mac ${WD}/MAC_bins2/MAC_bins_${nsim}.txt \
#    -o ${WD}/MAC_bins2/${pop}/Expected_MAC_${nsim}_${pop}_100.txt \
#    -N $((2*$ncase)) \
#    --nvar_target_data ${WD}/MAC_bins2/${pop}/chr19_block37_${pop}_nvar_target_data.txt \
#    --afs_target_data ${WD}/MAC_bins2/${pop}/chr19_block37_${pop}_AFS_target_data.txt \
#    --reg_size 19.029


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

# PRUNE HAPLOTYPES FOR BOTH SCENARIOS

# loop through the indices of the pcase percentages
for p in "${!pcase[@]}"; do
if [[ $p -eq 0 ]]; then

echo "Pruning down to ${pcase[$p]}% functional and 100% synonymous variants..."

# prune functional and synonymous variants down to pcase% functional / 100% synonymous and remove the rows of zeros (use the z flag)
python3 /home/math/murphjes/raresim2/sim.py \
    -m ${WD}/Hapgen$((nsim/1000))K/${prefix}.controls.haps.sm \
    --functional_bins ${WD}/MAC_bins2/${pop}/Expected_MAC_${nsim}_${pop}_${pcase[$p]}_fun.txt \
    --synonymous_bins ${WD}/MAC_bins2/${pop}/Expected_MAC_${nsim}_${pop}_100_syn.txt \
    --stop_threshold 10 \
    -l ${WD}/Hapgen$((nsim/1000))K/${prefix}.copy.legend \
    -L ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.${pcase[$p]}fun.100syn.legend \
    -H ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.${pcase[$p]}fun.100syn.haps.gz \
    -z

# convert the entire sample to a sparse matrix
python3 /home/math/murphjes/raresim2/convert.py \
    -i ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.${pcase[$p]}fun.100syn.haps.gz \
    -o ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.${pcase[$p]}fun.100syn.haps.sm

else

echo "Pruning down to ${pcase[$p]}% functional variants..."

# prune functional variants down to pcase% in the entire sample 
# (the legend file is not output when using the z flag, only the pruned variants file)
python3 /home/math/murphjes/raresim2/sim.py \
    -m ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.${pcase[$((p-1))]}fun.100syn.haps.sm \
    --f_only ${WD}/MAC_bins2/${pop}/Expected_MAC_${nsim}_${pop}_${pcase[$p]}_fun.txt \
    --stop_threshold 10 \
    -l ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.${pcase[0]}fun.100syn.legend \
    -L ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.${pcase[$p]}fun.100syn.legend \
    -H ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.${pcase[$p]}fun.100syn.haps.gz

# convert the entire sample to a sparse matrix
python3 /home/math/murphjes/raresim2/convert.py \
    -i ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.${pcase[$p]}fun.100syn.haps.gz \
    -o ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.${pcase[$p]}fun.100syn.haps.sm

fi

done

# EXTRACT CASES FOR SAME DIRECTION OF EFFECT SCENARIO

# loop through the indices of the percentages for the same direction of effect scenario
for s in "${!same[@]}"; do

# extract the power cases for the same direction of effect scenario
python3 /home/math/murphjes/raresim2/extract.py \
    -i ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.${same[$s]}fun.100syn.haps.gz \
    -o ${WD}/Cases2/${prefix}.${ncase}.power.cases.same.${same[$s]}fun.100syn.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

# remove unnecessary files
rm ${WD}/Cases2/${prefix}.${ncase}.power.cases.same.${same[$s]}fun.100syn.haps-remainder.gz

done

# PRUNE / SUBSET FOR OPPOSITE DIRECTIONS OF EFFECT SCENARIO

# loop through the indices of the percentages for the opposite directions of effect scenario
for o in "${!opp[@]}"; do

echo "Pruning down from ${opp[$o]}% to 100% functional variants..."

# prune functional variants down to 100% in the entire sample
python3 /home/math/murphjes/raresim2/sim.py \
    -m ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.${opp[$o]}fun.100syn.haps.sm \
    --f_only ${WD}/MAC_bins2/${pop}/Expected_MAC_${nsim}_${pop}_100_fun.txt \
    --stop_threshold 10 \
    -l ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.${pcase[0]}fun.100syn.legend \
    -L ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.100fun.100syn.${opp[$o]}.legend \
    -H ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.100fun.100syn.${opp[$o]}.haps.gz

# record the pruned variants (make sure to account for the header or lack thereof)
awk -F'\t' -v OFS='\t' 'NR==FNR {if (FNR > 1) ids[$1]=1; next} FNR==1 {print $0, "protected"; next} {print $0, ($1 in ids) ? 1 : 0}' \
${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.100fun.100syn.${opp[$o]}.legend-pruned-variants \
${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.${pcase[0]}fun.100syn.legend > ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.${opp[$o]}fun.100syn.protected.legend 

# extract the type I error cases (and power cases for the opposite direction of effect scenario)
python3 /home/math/murphjes/raresim2/extract.py \
    -i ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.100fun.100syn.${opp[$o]}.haps.gz \
    -o ${WD}/Cases2/${prefix}.${ncase}.t1e.cases.100fun.100syn.${opp[$o]}.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

# move the controls for the type I error scenario
mv ${WD}/Cases2/${prefix}.${ncase}.t1e.cases.100fun.100syn.${opp[$o]}.haps-remainder.gz ${WD}/Controls2/${prefix}.${ncase}.controls.same.100fun.100syn.${opp[$o]}.haps-remainder.gz

echo "Pruning down from ${opp[$o]}% to 100% functional variants for the opposite direction of effect scenario..."

# prune functional variants down to 100% in the entire sample again excluding the previously pruned variants from the entire sample
python3 /home/math/murphjes/raresim2/sim.py \
    -m ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.${opp[$o]}fun.100syn.haps.sm \
    --f_only ${WD}/MAC_bins2/${pop}/Expected_MAC_${nsim}_${pop}_100_fun.txt \
    --stop_threshold 10 \
    -l ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.${opp[$o]}fun.100syn.protected.legend \
    -L ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.opp.100fun.100syn.${opp[$o]}.legend \
    -H ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.opp.100fun.100syn.${opp[$o]}.haps.gz \
    --keep_protected

# extract the controls for the opposite direction of effect scenario
python3 /home/math/murphjes/raresim2/extract.py \
    -i ${WD}/Hapgen$((nsim/1000))K_pruned2/${prefix}.${nsim}.all.opp.100fun.100syn.${opp[$o]}.haps.gz \
    -o ${WD}/Controls2/${prefix}.${ncase}.controls.opp.100fun.100syn.${opp[$o]}.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

# remove unnecessary files
rm ${WD}/Controls2/${prefix}.${ncase}.controls.opp.100fun.100syn.${opp[$o]}.haps-sample.gz

done

done

# power cases (same): chr19.block37.${pop}.sim${rep}.${ncase}.power.cases.same.${pcase[$p]}fun.100syn.haps-sample.gz
# type I error cases & power cases (opp): chr19.block37.${pop}.sim${rep}.${ncase}.t1e.cases.100fun.100syn.${opp[$o]}.haps-sample.gz
# controls (same): chr19.block37.${pop}.sim${rep}.${ncase}.t1e.cases.100fun.100syn.${opp[$o]}.haps-remainder.gz
# controls (opp): chr19.block37.${pop}.sim${rep}.${ncase}.controls.opp.100fun.100syn.${opp[$o]}.haps-remainder.gz
