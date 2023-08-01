#!/bin/bash

while getopts p:o:c:s:l:i: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        p) partis=${OPTARG};;
        o) output=${OPTARG};;
        c) clones=${OPTARG};;
        s) shm=${OPTARG};;
        l) leaves=${OPTARG};;
        i) sim=${OPTARG};;
    esac
done

##### Step 2: Simulate Sequences #####
#test
## simulate and rearrange from scratch:


$partis simulate --simulate-from-scratch --n-sim-events ${clones} --scratch-mute-freq ${shm//_/.} --n-leaves ${leaves} --outfname $output/clones_${clones}_shm_${shm}_leaves_${leaves}_sim_${sim}.yaml --debug 1


#clonal_families=${output##*sim_}
#clonal_families=${clonal_families%/*}
#shm=${output##*/}
#echo ${shm}
#echo ${clonal_families}

#set x for number of desired simulations $(seq 1 x)
#for i in 0.01 0.05 0.1 0.2 0.3; do
 # echo "$i";
  #$partis simulate --simulate-from-scratch --n-sim-events ${name} --scratch-mute-freq ${i} --outfname $output/sim_${i//./_}.yaml --debug 1
#done


#simulate from existing
#$partis simulate --parameter-dir $output --n-sim-events ${i} --outfname $output/simulations/sim_${i}.yaml --debug 1
#for example file add
#--min-observations-per-gene 5 --n-leaf-distribution geometric;


