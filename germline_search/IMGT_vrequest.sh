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

for fasta in $directory/*.fasta; do
  name=${fasta%.fasta}
  name=${name##*/}
  ${vquest} --species human --receptorOrLocusType IG --fileSequences ${fasta} -o ${output} && mv ${output}/001 ${output}/${name}
done