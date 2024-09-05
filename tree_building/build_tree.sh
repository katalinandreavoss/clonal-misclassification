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

${raxml} -model GTR -msa ${directory}/clean_aligned.fasta -seed 42 -prefix ${output}/mega_tree_ --search ML tree search --blopt nr_safe --threads auto{64}
#for fasta in $directory/*_aligned.fasta; do
#  name=${fasta%_aligned.fasta}
#  name=${name##*/}
#  lines=$(cat ${fasta} | grep ">" | wc -l)
#  if [ ${lines} -gt 4 ]; then
#    ${raxml} -model GTR -msa ${fasta} -seed 42 -prefix ${output}/${name}_tree_ --search ML tree search
#  fi
#done


