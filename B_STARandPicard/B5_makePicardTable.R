#!/usr/bin/env Rscript
#$ -cwd
#$ -o  makepicardTab.$JOB_ID.out
#$ -j y
#  Resources requested:
#$ -l h_data=8G,h_rt=01:00:00
#$ -V
#  Email address to notify
#$ -M videlock@mail
#$ -m bea

options(stringsAsFactors=FALSE)

file.patterns=c("*.metricsAlign$","*rna_metrics$","*.metricsDup$")
file.path="picardout"

picList<-lapply(file.patterns, function(x){
  ff <- list.files( path = file.path , pattern = x, full.names = TRUE )
  pic.files <- lapply( ff, read.table, skip = 6,sep = "\t",quote="",header=TRUE,nrows=1)
  pic.df<-as.data.frame( sapply( pic.files, function(x) x[ 1, ,drop=F] ) )
  names(pic.df)<-gsub("picardout/","", gsub("\\..+","",ff))
  pic.df<-as.data.frame(t(pic.df))
  return(pic.df)
})
  
picard<-do.call(cbind,picList)

ctinfo<-readRDS("star_countinfo.rds")
ctinfo<-as.data.frame(t(ctinfo))
ctinfo<-ctinfo[rownames(picard),]
picard<-cbind(ctinfo,picard)


picard<-picard[,c(
  "N_unmapped","N_multimapping","N_noFeature","N_ambiguous",
  #AlignmentSummaryMetrics
  "PF_READS", "PF_HQ_ALIGNED_READS",
  "PF_INDEL_RATE",  "PF_HQ_ERROR_RATE", 
  "PF_ALIGNED_BASES", "PF_HQ_ALIGNED_BASES", "PF_HQ_ALIGNED_Q20_BASES", 
  "PF_MISMATCH_RATE", 
  
  #rnaseqmetrics
  "PF_BASES", "CODING_BASES", "PCT_CODING_BASES", 
  "UTR_BASES", "PCT_UTR_BASES",
  "INTRONIC_BASES", "PCT_INTRONIC_BASES",
  "INTERGENIC_BASES", "PCT_INTERGENIC_BASES", 
  "CORRECT_STRAND_READS", "PCT_CORRECT_STRAND_READS",
  "INCORRECT_STRAND_READS",
  "PCT_MRNA_BASES", "PCT_USABLE_BASES",
  "MEDIAN_CV_COVERAGE", "MEDIAN_3PRIME_BIAS",  
  
  #DuplicationMetrics
  "UNPAIRED_READ_DUPLICATES",  "PERCENT_DUPLICATION"
)]

saveRDS(picard,file = "picardTable.rds")

