#!/bin/bash

while getopts p:s:o: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        p) partis=${OPTARG};;
        s) simu=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

#clonal_families=${output##*sim_}
#clonal_families=${clonal_families%/*}
#shm=${output##*/}

python analyze_partis_output/yaml_to_families_new.py ${simu} ${partis} ${output}/

#for simu in $directory/sim_*.yaml; do
#  name=${simu%.yaml}
#  name=${name##*/}
#  echo ${name}
#  python analyze_partis_output/yaml_to_families_new.py ${simu} ${partis} ${output}/
#done