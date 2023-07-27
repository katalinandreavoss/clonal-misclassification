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
  lines=$(cat ${fasta} | grep ">" | wc -l)
  if [ ${name} != "all" ] && [ ${lines} -gt 4 ] ; then
    ${raxml} --ancestral -model GTR -msa ${fasta} --tree ${directory}/tree_files/${name}_tree_.raxml.bestTree -prefix ${output}/${name}
  fi
done

for i in ${output}/*.raxml.ancestralStates; do  sed '/^N/ s/./>&/' $i > $i.txt; done
for i in ${output}/*.raxml.ancestralStates.txt; do  sed '/>/h; //,$ { //!H; }; $!d; x' $i > $i.txt; done
cat ${output}/*txt.txt > ${output}/root_naive.txt #all together
ls ${output}/*txt.txt > ${output}/list.txt #file list

