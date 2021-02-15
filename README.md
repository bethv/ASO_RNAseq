# ASO_RNAseq

## ASO RNA-seq Pipeline

Code exceuted in the following order:

### A. Trim FASTQs
  1. bbduk - [A1_bbdukQS.sh](A_TrimFastqs/A1_bbdukQS.sh)
  2. FastQC - [A2_runFastQC.sh](A_TrimFastqs/A2_runFastQC.sh)

### B. Align reads with STAR and run Picard tools
  1. generate index - [B1_genomeGenMM22](B_STARandPicard/B1_genomeGenMM22.sh)
  2. align - [B2_runStar.sh](B_STARandPicard/B2_runStar.sh)
  3. samtools - [B3_samtoolsAndCountTable.sh](B_STARandPicard/B3_samtoolsAndCountTable.sh)
  4. count table - [makeCountTable.R](B_STARandPicard/makeCountTable.R)
  5. Picard tools - [B4_runPicard.sh](B_STARandPicard/B4_runPicard.sh)
  6. Picard table - [B5_makePicardTable.R](B_STARandPicard/B5_makePicardTable.R)

### C. QC and Normalization
  * outlier detection - [outlierDetection.R](C_QCandNormalization/outlierDetection.R)
  * 
  * colon, both 1 and 3 months - [dc_1and3m.R](C_QCandNormalization/dc_1and3m.R)
  * colon, 1 month - [dc_1m.R](C_QCandNormalization/dc_1m.R)
  * colon, 3 months - [dc_3m.R](C_QCandNormalization/dc_3m.R)
  * striatum, both 1 and 3 months - [str_1and3m.R](C_QCandNormalization/str_1and3m.R)
  * striatum, 1 month - [str_1m.R](C_QCandNormalization/str_1m.R)
  * striatum, 3 months - [str_3m.R](C_QCandNormalization/str_3m.R)
