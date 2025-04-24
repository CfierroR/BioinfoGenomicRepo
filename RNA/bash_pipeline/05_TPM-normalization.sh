#!/bin/bash

input_dir=$1
output_dir=$2

#for input_file in ${input_dir} 
find ${input_dir} -name "*.counts" | while read input_file;
do
	echo "Normalize ${input_file} file"
	outname=$(basename $input_file | sed 's/\.counts//g')
	SF=$(awk 'BEGIN{FS=OFS="\t";sum=0}{if(NR>2){sum=sum+$7/($6/1000)}}END{print sum/1000000}' ${input_file})
	echo "Scaling Factor for ${outname} : ${SF}"
	awk -v SF=${SF} 'BEGIN{FS=OFS="\t"}{if(NR>2){print $1,$6,$7,($7/($6/1000))/SF}}' ${input_file} > ${output_dir}/${outname}".tpm.txt"

done
