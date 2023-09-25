#!/bin/bash

while getopts p:d: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        p) ptp=${OPTARG};;
        d) directory=${OPTARG};;
    esac
done

for dir in ${directory}*/; do
    echo "$dir"
    python ${ptp} -t ${dir}tree_files/mega_tree_.raxml.bestTree -o ${dir}mega
done