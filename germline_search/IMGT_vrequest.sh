#!/bin/bash

while getopts d:o:v: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
        v) vquest=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

for fasta in $directory/*_aligned.fasta; do
  name=${fasta%_aligned.fasta}
  name=${name##*/}
  cd ${vquest}
  python -m vquest --species human --receptorOrLocusType IG --fileSequences ${fasta} -o ${output} && mv ${output}/001 ${output}/${name}
done