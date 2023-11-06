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
#if ancestral sequences for correct: add tree_files to raxml tree

rm -f ${output}/"*_.raxml.ancestralStates.txt"
rm -f ${output}/"*_.raxml.ancestralStates.txt.txt"
rm -f ${output}/"root_naive.txt"
rm -f ${output}/"*[^rerooted].raxml.ancestralStates.txt.txt"
rm -f ${output}/"*[^rerooted].raxml.ancestralStates.txt"
rm -f ${output}/*first*
rm -f ${output}/"*.raxml.ancestralStates.txt"
rm -f ${output}/"*.raxml.ancestralStates.txt.txt"

for fasta in $directory/family*_aligned.fasta; do
  name=${fasta%_aligned.fasta}
  name=${name##*/}
  lines=$(cat ${fasta} | grep ">" | wc -l)
  if [ ${name} != "clean" ] && [ ${lines} -gt 4 ] ; then
    ${raxml} --ancestral -model GTR -msa ${fasta} --tree ${output}/${name}_tree_.raxml.bestTree -prefix ${output}/${name}
  fi
done

for i in ${output}*[^rerooted].raxml.ancestralStates; do  sed '/^N/ s/./>&/' $i > $i.txt; done
for i in ${output}*[^rerooted].raxml.ancestralStates.txt; do  sed '/>/h; //,$ { //!H; }; $!d; x' $i > $i.txt; done

for fasta in ${output}/*[^rerooted].raxml.ancestralStates.txt.txt; do
  name=${fasta%.raxml.ancestralStates.txt.txt}
  name=${name##*/}
  echo -e -n ${name}"\t" >> ${output}/root_naive.txt
  awk '{print $2}' ${fasta} >> ${output}/root_naive.txt #all together
done

ls ${output}/*[^rerooted].raxml.ancestralStates.txt.txt > ${output}/list.txt #file list

