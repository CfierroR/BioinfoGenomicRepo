#!/bin/bash

# Script to run fastQC on fastq files.
# Input: a folder with the fastq files

input_dir=$1
echo "Running fastqc on the files: $(ls ${input_dir})"

out_dir="../outputs/01_fastqc/"
mkdir -p ${out_dir}
# Create tmp folder
mkdir -p ${out_dir}/tmp_dir

fastqc -o ${out_dir} --format fastq --threads 30 --dir ${out_dir}/tmp_dir ${input_dir}/*.fastq.gz

# delete tmp dir
rm -rf ${out_dir}/tmp_dir
