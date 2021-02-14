#!/bin/bash
#$ -cwd
#$ -j y
#  Resources requested:
#$ -l h_data=8G,h_rt=04:00:00
#$ -V


function runFastQC(){
	../FastQC/fastqc *.fastq.gz --outdir=../fastqcOut

}


tar -xvf trimmed_clean_fastqs.tar.gz 

cd trimmed_clean && runFastQC