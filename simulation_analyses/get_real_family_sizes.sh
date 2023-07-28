#!/bin/bash

while getopts d: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
    esac
done

clonal_families=${directory##*sim_}
clonal_families=${clonal_families%/*}
shm=${directory##*/}

for fasta in $directory/family*.fasta; do
  name=${fasta%.fasta}
  name=${name##*/}
  if [[ ${name} != *"_aligned" ]] ; then
    echo $shm
    echo ${shm} "," $(cat ${fasta} | wc -l | awk '{print $1/2}')| tr -d '\n' >> $directory/family_sizes.txt
    echo "" >> $directory/family_sizes.txt
  fi
done
