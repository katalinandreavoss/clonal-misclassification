#!/bin/bash

while getopts d:o: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

for fasta in $directory*.fasta; do
  name=${fasta%.fasta}
  name=${name##*/}
  echo $name
  clustalo -i $fasta -t DNA -o ${output}/${name}_aligned.fasta
done