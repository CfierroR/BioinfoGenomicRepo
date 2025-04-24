#!/bin/bash

input_dir=$1
output_dir=$2
s_opt=$3

function ftcount {
    bam=$1
    gtf=$2
    out=$3
    s_option=$4
    outname=$(basename $bam | sed 's/Aligned.sortedByCoord.out.bam//g')
    featureCounts -p -T 10 -s $s_option -a ${gtf} -o "${out}/opt_${s_opt}/${outname}.counts" ${bam}
}

gtf_mm10="/home/resources/genomes/Genomes/ensembl/hg38/Homo_sapiens.GRCh38.104.chr.gtf"

find ${input_dir} -name "*.sortedByCoord.out.bam" | while read file;
do
    #outdir=$(echo $file | cut -f1-4 -d'/')
    mkdir -p ${output_dir}
    #if [[ -n $(echo $file | grep "N1"|"N2" ) ]]; then
    #    s_opt=0
    #else
    #    s_opt=0
    #fi
    mkdir -p "${output_dir}/opt_${s_opt}"
    ftcount ${file} ${gtf_mm10} ${output_dir} ${s_opt}


done
