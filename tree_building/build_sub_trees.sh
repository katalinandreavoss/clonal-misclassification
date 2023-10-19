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

for fasta in $directory/family_*_aligned.fasta; do
  name=${fasta%_aligned.fasta}
  name=${name##*/}
  lines=$(cat ${fasta} | grep ">" | wc -l)
  if [ ${lines} -gt 4 ]; then
    ${raxml} -model GTR -msa ${fasta} -seed 42 -prefix ${output}/${name}_tree_ --search ML tree search
  fi
done
echo "build_tree done" > ${output}/build_tree.txt