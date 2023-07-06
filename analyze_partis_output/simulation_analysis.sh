#!/bin/bash

while getopts p:d:o: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        p) partis=${OPTARG};;
        d) directory=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

for simu in $directory/*sim_*.yaml; do
  name=${simu%.yaml}
  name=${name##*/}
  python analyze_partis_output/yaml_to_families_new.py ${simu} ${partis} ${output}/${name}
done