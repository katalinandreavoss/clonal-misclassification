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
  cd ${vquest}
  python -m vquest --species human --receptorOrLocusType IG --fileSequences ${fasta} -o ${output} && mv ${output}/001 ${output}/${name}
  [ -d ${output}/${name}/001 ] && cd ${output}/${name}/001 && mv * ${output}/${name} && rm -rf ${output}/${name}/001 #for some reason the last simulation gets stores incorrectly, this fixes it
  echo "one loop done";
done

echo "completely done"
