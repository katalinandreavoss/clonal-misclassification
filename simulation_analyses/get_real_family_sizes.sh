#!/bin/bash

while getopts d:o:c:s:l:i: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
        o) output=${OPTARG};;
        c) clones=${OPTARG};;
        s) shm=${OPTARG};;
        l) leaves=${OPTARG};;
        i) sim=${OPTARG};;
    esac
done

for fasta in $(ls $directory| grep -E "family_[0-9]+.fasta"); do
  lines=$(cat $directory/${fasta} | grep ">" | wc -l)
  if [ ${lines} -gt 1 ]; then
    echo ${clones} "," ${shm} "," ${leaves} "," ${sim}  "," ${lines} | tr -d '\n' >> $directory/family_sizes.txt
    echo "" >> $directory/family_sizes.txt
  fi
done
