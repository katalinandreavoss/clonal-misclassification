#!/bin/bash

while getopts p:f:o: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        p) partis=${OPTARG};;
        f) fasta=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

# partition the fasta file used in cache-parameters to create multi-hmm folder in the parameter directory
$partis partition --infname $fasta --outfname $output/pd.yaml --parameter-dir $output --count-parameters --min-observations-per-gene 5

##### Step 2: Simulate Sequences #####

## simulate and rearrange from scratch:

#set x for number of desired simulations $(seq 1 x)
for i in $(seq 1 5); do
  echo "$i";
  $partis simulate --parameter-dir $output --n-sim-events ${i} --outfname sim_${i}.yaml;
done
