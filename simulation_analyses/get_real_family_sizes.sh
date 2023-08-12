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

for fasta in $directory/family*.fasta; do
  echo ${clones} "," ${shm} "," ${leaves} "," ${sim}  "," $(cat ${fasta} | wc -l | awk '{print $1/2}')| tr -d '\n' >> $directory/family_sizes.txt
  echo "" >> $directory/family_sizes.txt
done
