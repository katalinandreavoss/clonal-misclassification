#!/bin/bash

while getopts r:o:f: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        r) raxml=${OPTARG};;
        o) output=${OPTARG};;
        f) filename=${OPTARG};;
    esac
done

${raxml} -model GTR -msa ${filename} -seed 42 -prefix ${output}/mega_tree_ --search ML tree search --blopt nr_safe --threads auto{64} --redo --log debug
#for fasta in $directory/*_aligned.fasta; do
#  name=${fasta%_aligned.fasta}
#  name=${name##*/}
#  lines=$(cat ${fasta} | grep ">" | wc -l)
#  if [ ${lines} -gt 4 ]; then
#    ${raxml} -model GTR -msa ${fasta} -seed 42 -prefix ${output}/${name}_tree_ --search ML tree search
#  fi
#done


