#!/bin/bash

while getopts d:f: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
        f) filename=${OPTARG};;
    esac
done

clustalo -i ${filename}.fasta -t DNA --outfmt=fasta -o ${filename}_aligned.fasta -v
