#!/bin/bash

while getopts d: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
    esac
done

clustalo -i ${directory}/clean.fasta -t DNA --outfmt=fasta -o ${directory}/clean_aligned.fasta
