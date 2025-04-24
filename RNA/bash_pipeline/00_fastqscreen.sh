#!/bin/bash

# Script to run fastQC on fastq files.
# Input: a folder with the fastq files

input_dir=$1
echo "Running fastqc on the files: $(ls ${input_dir})"

out_dir="../outputs/00_fastqScreen/"
mkdir -p ${out_dir}
# Create tmp folder
mkdir -p ${out_dir}/tmp_dir

/home/resources/programs/fastq_screen/FastQ-Screen-0.14.1/fastq_screen -o ${out_dir} --threads 30 ${input_dir}/*.fastq.gz

# delete tmp dir
rm -rf ${out_dir}/tmp_dir
