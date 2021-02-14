#!/bin/bash
#$ -l h_rt=5:00:00,h_data=32G
#$ -cwd -V
#$ -j y
#$ -o poststar$JOB_ID.out
#$ -M videlock@mail
#$ -m bea

. /u/local/Modules/default/init/modules.sh
module load STAR/2.7.0e

samples=$(cat samps)
nthreads=8
mkdir -p readcounts
mkdir -p bamfiles


cp starout/*ReadsPerGene.out.tab readcounts

cp starout/*Aligned.sortedByCoord.out.bam bamfiles

for i in ${samples}; do
	STAR --runThreadN ${nthreads} --runMode inputAlignmentsFromBAM \
	--inputBAMfile bamfiles/${i}_Aligned.sortedByCoord.out.bam \
	--outWigType bedGraph --outWigStrand Unstranded
done


module load samtools
cd bamfiles
for bamfile in *Aligned.sortedByCoord.out.bam ; do samtools index ${bamfile}; done
cd ../


module load R/3.5.1

R CMD BATCH makeCountTable.R