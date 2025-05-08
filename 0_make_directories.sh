#!/bin/bash

Nsim=10000

mkdir ./results
mkdir ./datasets
mkdir ./datasets/Hapgen$((Nsim/1000))K
mkdir ./datasets/Hapgen$((Nsim/1000))K_pruned
mkdir ./datasets/Cases
mkdir ./datasets/Controls

pops=(AFR EAS NFE SAS)

for pop in "${pops[@]}"; do

mkdir ./results/${pop}
mkdir ./datasets/Hapgen$((Nsim/1000))K/${pop}
mkdir ./datasets/Hapgen$((Nsim/1000))K_pruned/${pop}
mkdir ./datasets/Cases/${pop}
mkdir ./datasets/Controls/${pop}

for n in {1..10}; do

mkdir ./datasets/Hapgen$((Nsim/1000))K/${pop}/Round${n}
mkdir ./datasets/Hapgen$((Nsim/1000))K_pruned/${pop}/Round${n}
mkdir ./datasets/Cases/${pop}/Round${n}
mkdir ./datasets/Controls/${pop}/Round${n}

done

done


