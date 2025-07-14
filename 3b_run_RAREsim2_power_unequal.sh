#!/bin/bash

#singularity shell /storage/singularity/mixtures.sif

nsim=10000 #number of individuals
num=19 #chromosome number
b=37 #block number
ncase=5000 #number of cases (also number of controls)
same=(160 140 120) #pcase necessary for same direction of effect scenario

# 75%/25% split of positive/negative directions of effect
opp1=(145 130 115) #more causal variants in the cases
opp2=(115 110 105)

#WD=/data001/projects/murphjes/RAREsim2/datasets

dir=/home/math/murphjes/raresim2 #RAREsim2 directory

# CREATE EXPECTED MAC BINS (need target data for all populations)

pcase=("${opp1[@]}" "${opp2[@]}")

for p in "${!pcase[@]}"; do

python3 ${dir}/expected_vars.py \
    --mac ./input/mac_bins/MAC_bins_${nsim}.txt \
    -o ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_${pcase[$p]}.txt \
    -N $((2*$ncase)) \
    --nvar_target_data ./input/mac_bins/${pop}/chr19_block37_${pop}_nvar_target_data.txt \
    --afs_target_data ./input/mac_bins/${pop}/chr19_block37_${pop}_AFS_target_data.txt \
    --w_fun $(awk "BEGIN { printf \"%.1f\", ${pcase[$p]} / 100 }") \
    --reg_size 19.029

rm ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_${pcase[$p]}_syn.txt

done

start=$(($end-99))

# loop through the simulation replicates
for i in $(eval echo "{$start..$end}")
do 

echo "Simulation Replicate: ${i}"

# determine n from i
n=$(( (i - 1) / 1000 + 1 ))

prefix=${pop}/Round${n}/chr${num}.block${b}.${pop}.sim${i}

# PRUNE / SUBSET HAPLOTYPES FOR THE OPPOSITE DIRECTIONS OF EFFECT SCENARIO

# loop through the indices of the percentages for the opposite directions of effect scenario
for o in "${!opp1[@]}"; do

echo "Pruning down from ${same[$o]}% to ${opp1[$o]}% functional variants..."

# prune functional variants down from same% to opp1% in the entire sample
python3 ${dir}/sim.py \
    -m ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.same.${same[$o]}fun.100syn.haps.sm \
    --f_only ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_${opp1[$o]}_fun.txt \
    --stop_threshold 10 \
    -l ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${same[0]}fun.100syn.legend \
    -L ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${opp1[$o]}fun2.100syn.legend \
    -H ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.opp.${opp1[$o]}fun2.100syn.haps.gz

# convert the entire sample to a sparse matrix
python3 ${dir}/convert.py \
    -i ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.opp.${opp1[$o]}fun2.100syn.haps.gz \
    -o ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.opp.${opp1[$o]}fun2.100syn.haps.sm

# record the pruned variants (make sure to account for the header or lack thereof)
awk -F'\t' -v OFS='\t' 'NR==FNR {if (FNR > 1) ids[$1]=1; next} FNR==1 {print $0, "protected"; next} {print $0, ($1 in ids) ? 1 : 0}' \
./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${opp1[$o]}fun2.100syn.legend-pruned-variants \
./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${same[0]}fun.100syn.legend > ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${opp1[$o]}fun2.100syn.protected.legend 

# extract the power cases for the opposite directions of effect scenario
python3 ${dir}/extract.py \
    -i ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.opp.${opp1[$o]}fun2.100syn.haps.gz \
    -o ./datasets/Cases/${prefix}.${ncase}.power.cases.opp.${opp1[$o]}fun2.100syn.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

echo "Pruning down from ${same[$o]}% to ${opp2[$o]}% functional variants, excluding protected variants..."

# prune functional variants down from same% to opp2% in the entire sample again excluding the previously pruned variants
python3 ${dir}/sim.py \
    -m ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.same.${same[$o]}fun.100syn.haps.sm \
    --f_only ./input/mac_bins/${pop}/Expected_MAC_${nsim}_${pop}_${opp2[$o]}_fun.txt \
    --stop_threshold 10 \
    -l ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.${opp1[$o]}fun2.100syn.protected.legend \
    -L ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.opp.${opp2[$o]}fun2.100syn.legend \
    -H ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.opp.protected.${opp2[$o]}fun2.100syn.haps.gz \
    --keep_protected

# extract the controls for the opposite direction of effect scenario
python3 ${dir}/extract.py \
    -i ./datasets/Hapgen$((nsim/1000))K_pruned/${prefix}.${nsim}.all.opp.protected.${opp2[$o]}fun2.100syn.haps.gz \
    -o ./datasets/Controls/${prefix}.${ncase}.controls.opp.${opp2[$o]}fun2.100syn.haps.gz \
    -n $((2*$ncase)) \
    --seed $i

# remove unnecessary files
rm ./datasets/Controls/${prefix}.${ncase}.controls.opp.${opp2[$o]}fun2.100syn.haps-sample.gz

done

done

# power cases (opp): chr19.block37.${pop}.sim${rep}.${ncase}.power.cases.opp.${opp1[$o]}fun2.100syn.haps-sample.gz
# controls (opp): chr19.block37.${pop}.sim${rep}.${ncase}.controls.opp.${opp2[$o]}fun2.100syn.haps-remainder.gz
