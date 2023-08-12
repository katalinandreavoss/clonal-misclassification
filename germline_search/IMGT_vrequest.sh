#!/bin/bash

while getopts f:o:v: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        f) fasta=${OPTARG};;
        v) vquest=${OPTARG};;
        o) output=${OPTARG};;
    esac
done


cd ${vquest}
python -m vquest --species human --receptorOrLocusType IG --fileSequences ${fasta} -o ${output}
