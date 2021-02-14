#!/bin/bash
#$ -l h_rt=8:00:00,h_data=32G,exclusive
#$ -cwd -V
#$ -N mm22stgi
#$ -j y
#$ -o gi.$JOB_ID.out
#$ -M videlock@mail
#$ -m bea

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz -O mm22.gtf.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/GRCm38.primary_assembly.genome.fa.gz -O mm22.fa.gz

gunzip mm22.fa.gz
gunzip mm22.gtf.gz

nthreads=8

mkdir -p genomeIndex


. /u/local/Modules/default/init/modules.sh
module load STAR/2.7.0e

STAR --runThreadN ${nthreads} --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles mm22.fa --sjdbGTFfile mm22.gtf --sjdbOverhang 64

