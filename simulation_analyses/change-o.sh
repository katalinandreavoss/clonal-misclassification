#!/bin/bash

while getopts d:f:v: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
        f) fasta=${OPTARG};;
        v) vdj=${OPTARG};;
    esac
done

MakeDb.py imgt -i ${directory} -s ${fasta} --extended
CreateGermlines.py -d ${directory}_db-pass.tsv -g dmask -r ${vdj}/IGHV.fasta ${vdj}/IGHD.fasta ${vdj}/IGHJ.fasta 
DefineClones.py -d ${directory}_db-pass.tsv
