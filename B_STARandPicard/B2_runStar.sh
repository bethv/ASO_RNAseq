#!/bin/bash
#$ -l h_rt=5:00:00,h_data=32G,exclusive
#$ -cwd -V
#$ -t 1-7
#$ -j y
#$ -o star$JOB_ID.out
#$ -M videlock@mail
#$ -m bea

samples=$(cat samps${SGE_TASK_ID})
nthreads=8

mkdir -p starout


. /u/local/Modules/default/init/modules.sh
module load STAR/2.7.0e


for i in ${samples}; do
	STAR --runThreadN ${nthreads} --genomeDir genomeIndex \
	--readFilesIn trimmed_clean/${i}.fastq \
	--outFilterType BySJout --outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
	--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
	--outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx --quantMode GeneCounts \
	--outFileNamePrefix starout/${i}_
done

