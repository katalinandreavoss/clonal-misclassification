#!/bin/bash

while getopts p:o:c:s:l:i:b: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        p) partis=${OPTARG};;
        o) output=${OPTARG};;
        c) clones=${OPTARG};;
        s) shm=${OPTARG};;
        l) leaves=${OPTARG};;
        b) balance=${OPTARG};;
        i) sim=${OPTARG};;
    esac
done

##### Step 2: Simulate Sequences #####
#test
## simulate and rearrange from scratch:
if [[ ${balance} == "0_0" ]]; then
    $partis simulate --simulate-from-scratch --n-genes-per-region 6:2:1 --n-sim-events ${clones} --scratch-mute-freq ${shm//_/.} --n-leaves ${leaves} --outfname $output/simu.yaml --debug 1
else 
    $partis simulate --simulate-from-scratch --n-sim-events ${clones} --scratch-mute-freq ${shm//_/.} --n-leaves ${leaves} --root-mrca-weibull-parameter ${balance//_/.} --outfname $output/simu.yaml --debug 1
fi


#--only-genes IGHV1-18*01:IGHV1-18*03:IGHV1-2*01:IGHV1-2*02:IGHV1-24*01:IGHV1-3*01:IGHV1-46*03:IGHD1-26*01:IGHD2-15*01:IGHJ6*01

#clonal_families=${output##*sim_}
#clonal_families=${clonal_families%/*}
#shm=${output##*/}
#echo ${shm}
#echo ${clonal_families}

#set x for number of desired simulations $(seq 1 x)
#for i in 0.002 0.0025 0.003 0.004; do
 # echo "$i";
  #$partis simulate --simulate-from-scratch --n-sim-events ${name} --scratch-mute-freq ${i} --outfname $output/sim_${i//./_}.yaml --debug 1
#done


#simulate from existing
#$partis simulate --parameter-dir $output --n-sim-events ${i} --outfname $output/simulations/sim_${i}.yaml --debug 1
#for example file add
#--min-observations-per-gene 5 --n-leaf-distribution geometric;


