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

## simulate and rearrange from scratch:

#set x for number of desired simulations $(seq 1 x)
for i in $(seq 1 250); do
  echo "$i";
  $partis simulate --parameter-dir $output --n-sim-events 5 --outfname $output/sim_${i}.yaml --min-observations-per-gene 5 --n-leaf-distribution geometric --debug 1
done
#for example file add
#--min-observations-per-gene 5 --n-leaf-distribution geometric;


