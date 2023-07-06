#!/bin/bash

while getopts r:d:o: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        r) raxml=${OPTARG};;
        d) directory=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

for fasta in $directory/*_aligned.fasta; do
  name=${fasta%_aligned.fasta}
  name=${name##*/}
  ${raxml} -model GTR -msa ${fasta} -seed 42 -prefix ${output}_${name}_tree_ --search ML tree search
done


