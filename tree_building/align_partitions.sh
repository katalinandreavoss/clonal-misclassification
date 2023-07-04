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
  echo $fasta
  name=${fasta%.fasta}
  echo $name
  clustalo -i $fasta -t DNA -o ${name}_aligned.fasta
done