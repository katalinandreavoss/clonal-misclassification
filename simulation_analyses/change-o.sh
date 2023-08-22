#!/bin/bash

while getopts d:f: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
        f) fasta=${OPTARG};;
    esac
done

MakeDb.py imgt -i ${directory} -s ${fasta} --extended
DefineClones.py -d ${directory}_db-pass.tsv
