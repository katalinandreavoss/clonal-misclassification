#!/bin/bash

while getopts d:o:r: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
        r) RTG=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

for dir in $directory/*/; do
  file=${dir}3_Nt-sequences.txt
  name=${dir##*/}
  python ${RTG}/RevertToGermlineCmd.py ${file} ${RTG}/imgt_germlines.fasta "Homo sapiens" ${output}/${name}/germline.fasta ocf
done