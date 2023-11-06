#!/bin/bash

while getopts d: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
    esac
done

for fasta in $directory/clean_cln*.fastq.gz; do
    name=${fasta%.fastq.gz}
    seqtk seq -a ${fasta} >  ${name}.fasta
    new_name=$(head ${name}.fasta -n 1 | cut -c2-)
    mv ${name}.fasta ${directory}/${new_name}.fasta 
done

echo "extract fastas done" > ${directory}/mixcr_fastas.txt