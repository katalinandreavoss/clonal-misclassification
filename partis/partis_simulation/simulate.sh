#!/bin/bash

while getopts p:o: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        p) partis=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

##### Step 2: Simulate Sequences #####
#test
## simulate and rearrange from scratch:
clonal_families=${output##*sim_}
clonal_families=${clonal_families%/*}
shm=${output##*/}
echo ${shm}
echo ${clonal_families}

$partis simulate --simulate-from-scratch --n-sim-events ${clonal_families} --scratch-mute-freq ${shm//_/.} --n-leaves 50 --outfname $output/sim_${shm}.yaml --debug 1

#set x for number of desired simulations $(seq 1 x)
#for i in 0.01 0.05 0.1 0.2 0.3; do
 # echo "$i";
  #$partis simulate --simulate-from-scratch --n-sim-events ${name} --scratch-mute-freq ${i} --outfname $output/sim_${i//./_}.yaml --debug 1
#done


#simulate from existing
#$partis simulate --parameter-dir $output --n-sim-events ${i} --outfname $output/simulations/sim_${i}.yaml --debug 1
#for example file add
#--min-observations-per-gene 5 --n-leaf-distribution geometric;


