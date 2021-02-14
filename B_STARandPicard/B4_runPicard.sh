#!/bin/bash
#$ -cwd
#$ -o  log.$JOB_ID.out
#$ -j y
#  Resources requested:
#$ -l h_data=8G,h_rt=04:00:00
#$ -t 1-7
#  Email address to notify
#$ -M videlock@mail
#$ -m bea
#$ -V



samples=$(cat ../samps${SGE_TASK_ID})

. /u/local/Modules/default/init/modules.sh
module load picard_tools
module load R/3.5.1


## Create mm22 dictionary for Picard
# java -Xmx5G -jar $PICARD CreateSequenceDictionary R=../mm22.fa

for i in ${samples}; do
# 	java -Xmx4G -jar $PICARD CollectRnaSeqMetrics REF_FLAT=../mm22refFlat.txt STRAND=FIRST_READ_TRANSCRIPTION_STRAND I=${i}_Aligned.sortedByCoord.out.bam O=${i}.rna_metrics CHART=${i}rna_metrics.pdf
	#java -Xmx4G -jar $PICARD CollectAlignmentSummaryMetrics R=../mm22.fa I=${i}_Aligned.sortedByCoord.out.bam O=${i}.metricsAlign
	java -Xmx4G -jar $PICARD CollectGcBiasMetrics R=../mm22.fa I=${i}_Aligned.sortedByCoord.out.bam O=${i}.metricsGCFull S=${i}.metricsGCSum CHART=${i}.metricsGC.pdf
	#java -Xmx4G -jar $PICARD MarkDuplicates I=${i}_Aligned.sortedByCoord.out.bam O=${i}.markDup.bam M=${i}.metricsDup
done



