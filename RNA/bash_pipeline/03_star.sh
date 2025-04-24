#!/bin/bash

# This script executes STAR
# Usage:
#   ./02_star.sh [input_dir] [output_dir]
# where:
#   [input_dir] is a folder with the fastq files
#   [output_dir] is the destination folder

input_dir=$1
outdir=$2

function star {
	inputs=$1
	genomedir=$2
	gtf=$3 
	sampleid=$(basename $4)
	
	#outdir="../outputs/02_star"
	mkdir -p ${outdir}/${sampleid}

	echo "Analyzing samples: $sampleid"
        STAR --runMode alignReads --runThreadN 30 --genomeDir $genomedir --sjdbGTFfile $gtf --sjdbOverhang 149 --readFilesIn $inputs --twopassMode Basic \
	--outFileNamePrefix ${outdir}/${sampleid}/${sampleid} \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMattributes Standard \
	--quantMode GeneCounts
	 
 
}

ref_mm10="/home/resources/genomes/Genomes/STAR_hg38"
gtf_mm10="/home/resources/genomes/Genomes/ensembl/hg38/Homo_sapiens.GRCh38.104.chr.gtf"

find ${input_dir} -name "*.fq.gz" |  cut -f8 -d'/' | cut -f1,2 -d'_' | sort | uniq | while read id ;
do
	# If single end
	#if [[ $(find ${input_dir} -name "*.fq.gz" | grep $id | wc -l) == 1 ]];
	#then
	#	fastqs=$(find ${input_dir} -name "*.fq.gz" | grep $id)
	#else	# else paired-end
	fq1=$(find ${input_dir} -name "*.fq.gz" | grep $id | grep "R1")
	fq2=$(find ${input_dir} -name "*.fq.gz" | grep $id | grep "R2")
	fastqs=$(echo "$fq1 $fq2")
	echo "$fq1 $fq2"
	echo "==============================================="

#	fi

	echo "#FASTQs: ${fastqs}"
	echo "#ref genome: ${ref_mm10}"
	echo "#gtf genome: ${gtf_mm10}"
	echo "#id: ${id}"
	
	star "$fastqs" $ref_mm10 $gtf_mm10 $id
	echo "-----------------------------------------------"
	
done
