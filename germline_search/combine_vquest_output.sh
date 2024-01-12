#!/bin/bash

while getopts d: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
        d) directory=${OPTARG};;
    esac
done


for dir in $directory*/; do
    echo ${dir}
    tail -n+2 ${dir}/1_Summary.txt >> $directory/1_Summary.txt
    tail -n+2 ${dir}/2_IMGT-gapped-nt-sequences.txt >> $directory/2_IMGT-gapped-nt-sequences.txt
    tail -n+2 ${dir}/3_Nt-sequences.txt >> $directory/3_Nt-sequences.txt
    tail -n+2 ${dir}/6_Junction.txt >> $directory/6_Junction.txt
done

for file in $directory*.txt; do
    echo ${file}
    echo "$(cut --complement -f 1 ${file})" > ${file}
    echo "$(awk '{printf "%s\t%s\n",NR,$0}' ${file})" > ${file}
done



line_1=$(head -n 1 $directory/001/1_Summary.txt)
sed -i "1i${line_1}/" $directory/1_Summary.txt
#echo ${line_1} > $directory/1_Summary.txt

line_2=$(head -n 1 $directory/001/2_IMGT-gapped-nt-sequences.txt)
sed -i "1i${line_1}/" $directory/2_IMGT-gapped-nt-sequences.txt

line_3=$(head -n 1 $directory/001/3_Nt-sequences.txt)
sed -i "1i${line_1}/" $directory/3_Nt-sequences.txt

line_6=$(head -n 1 $directory/001/6_Junction.txt)
sed -i "1i${line_1}/" $directory/6_Junction.txt