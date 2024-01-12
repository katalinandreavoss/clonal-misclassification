#!/bin/bash

while getopts d: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
    esac
done
sleep 5
for fasta in $directory/*_combined.fasta; do
  name=${fasta%.fasta}
  name=${name##*/}
  clustalo -i $fasta -t DNA --outfmt=fasta -o ${directory}/${name}_aligned.fasta 
done