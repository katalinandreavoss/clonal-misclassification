#!/bin/bash

while getopts d: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
    esac
done
#change to and from fasta/tsv 
for fasta in $directory/family*.fasta; do
  name=${fasta%.tsv}
  name=${name##*/}
  lines=$(cat ${fasta} | wc -l)
  if [ ${name} != "naive" ] && [ ${lines} -gt 2 ]; then
    clustalo -i $fasta -t DNA --outfmt=fasta -o ${directory}/${name}_aligned.fasta
  fi
  if [ ${name} != "naive" ] && [ ${lines} -eq 2 ]; then
    cp $fasta ${directory}/${name}_aligned.fasta
  fi
done

echo "done" > $directory/align.txt
