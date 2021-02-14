rm(list = ls())
source("QCfunctions.R")
grp="dc1m"
outname<-paste0(grp,"f2")
# setup ---------
library(Biobase)
library(limma)
library(edgeR)
library(DESeq2)
library(matrixStats)
library(MASS)
library(pheatmap)
library(tidyverse)
library(gridExtra)
library(grid)

dir.create(outname,showWarnings = FALSE)
datdir<-paste0(outname,"/data/")
figdir<-paste0(outname,"/figs/")
dir.create(datdir,showWarnings = FALSE)
dir.create(figdir,showWarnings = FALSE)
dir.create(paste0(outname,"/tex/"),showWarnings = FALSE)

eset<-readRDS("data/rawCountEset.rds") 
eset$Reads<-as.numeric(gsub(",","",eset$Reads))

outliers=c("7648_dc","7675_dc")
eset<-eset[,!rownames(pData(eset))%in%outliers]
eset<-eset[,paste0(eset$tissue,eset$Time,"m")==grp]
dim(eset);table(eset$GT)

# filter/vsd ---------------
temp<-as.data.frame(exprs(eset))
temp$flag<-0
for (i in 1:ncol(exprs(eset))){
  temp$flag<-temp$flag+(temp[,i]>=5)
}
thirtyPercent=round(.3*dim(eset)[2],digits = 0)
thirtyPercent
table(temp$flag>=thirtyPercent) 
keep=temp$flag>=thirtyPercent
filteredEset<-eset[keep,]
dim(filteredEset)
filtEsetFile<-paste0(datdir,"filteredEset.rds")
saveRDS(filteredEset,filtEsetFile)

cds <- DESeqDataSetFromMatrix(countData = exprs(filteredEset),
                              colData = pData(filteredEset),
                              design = ~ GT)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
vsd <- getVarianceStabilizedData(cds)
normDat<-vsd
saveRDS(normDat,paste0(datdir,"NormCounts.rds"))
#make into eset
filtNormEset<-filteredEset
exprs(filtNormEset)<-normDat
filtNormEsetFile<-paste0(datdir,"filtNormEset.rds")
saveRDS(filtNormEset,filtNormEsetFile)

# density plots -------------------
pdf(paste0(figdir,"densityplots.pdf"),width = 7.5*2, height = 9*2)

lm=rbind(c(13,8,8,8),
         c(9,1,2,2),
         c(10,3,4,4),
         c(11,5,6,6),
         c(12,7,7,7))

layout(mat = lm,
       heights = c(1,6,6,6,8))

dlines(exprs(eset), main = paste0("n=",dim(eset)[1]))
dboxes(exprs(eset))
dlines(exprs(filteredEset),main = paste0("n=",dim(filteredEset)[1]))
dboxes(exprs(filteredEset))
dlines(exprs(filtNormEset),log = F,main = paste0("n=",dim(filtNormEset)[1]))
dboxes(exprs(filtNormEset),log = F)
plotDispEsts(cds)

par(mar = c(0,0,0,0))
plot.new()
text(0.5,0.5,paste0("striatum, n=",dim(eset)[2]),cex=3)

plot.new();text(0.5,0.5,"raw",cex=3)
plot.new();text(0.5,0.5,"filtered",cex=3)
plot.new();text(0.5,0.5,"VSD",cex=3)
dev.off()

pdf(paste0(figdir,"densityplots2.pdf"),width = 8, height = 9)
lm=rbind(c(7,1,2,2),
         c(8,3,4,4),
         c(9,5,6,6))
layout(lm)
dlines(exprs(eset), main = paste0("n=",dim(eset)[1]))
dboxes(exprs(eset))
dlines(exprs(filteredEset),main = paste0("n=",dim(filteredEset)[1]))
dboxes(exprs(filteredEset))
dlines(exprs(filtNormEset),log = F,main = paste0("n=",dim(filtNormEset)[1]))
dboxes(exprs(filtNormEset),log = F)

par(mar = c(0,0,0,0))
plot.new();text(0.5,0.5,"raw",cex=3)
plot.new();text(0.5,0.5,"filtered",cex=3)
plot.new();text(0.5,0.5,"VSD",cex=3)
dev.off()

pdf(paste0(figdir,"meanvar.pdf"),width = 5, height = 5)
plotDispEsts(cds)
dev.off()

# pre - ----------
conds<-pData(filtNormEset)
conds$sample<-rownames(conds)
conds$RIN<-as.numeric(conds$RIN)
conds<-conds[,c("sample","RIN","GT")]
conds2<-conds[,c("sample","GT")]




mdslist<-list()
mdslist$pre<-makemds.set(data = normDat)

mdsplot1v2(mdslist$pre,figdir,"MDS_pre.pdf")
mdsplotAll(mdslist$pre,figdir,"MDS_pre_all.pdf")


# picard --------------
picard<-readRDS("data/picardTable.rds")
picard<-picard[,c(
  "N_unmapped", "N_multimapping","N_noFeature","N_ambiguous", #star count info
  
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
  "MEDIAN_CV_COVERAGE", 
  
  #DuplicationMetrics
  "UNPAIRED_READ_DUPLICATES",  "PERCENT_DUPLICATION")]

picardgrp<-picard[rownames(conds),]
picardgrp<-as.data.frame(apply(picardgrp, 2, as.numeric))
rownames(picardgrp)<-rownames(conds);head(picardgrp)
target<-cbind(conds,picardgrp);dim(target)
n=dim(conds)[2]
target<-target[,!(colSums(is.na(target)) > 0)];dim(target)
temp1<-target[,-c(1:n)];print(dim(temp1))
temp1<-temp1[,colVars(as.matrix(temp1))!=0];dim(temp1)
temp1<-t(scale(temp1))
temp1<-as.data.frame(t(temp1))
target<-cbind(target[,1:n],temp1)

normExpr<-as.data.frame(exprs(filtNormEset))
unique(match(names(normExpr),target$sample)==seq(1,nrow(target),1))
tar1=target[,c((1+dim(conds)[2]):dim(target)[2])]
thisdat <- t(scale(tar1,scale=F))
PC.metadata <- prcomp(thisdat,center=F);
topPC1 <- PC.metadata$rotation[,1:10];
varexp <- (PC.metadata$sdev)^2 / sum(PC.metadata$sdev^2)
seqtopvar <- varexp[1:10]
colnames(topPC1) <- paste("Seq.PC\n",colnames(topPC1)," (",signif(100*seqtopvar[1:7],2),"%)",sep="")
target$Seq.PC1=as.numeric(topPC1[,1])
target$Seq.PC2=as.numeric(topPC1[,2])
target$Seq.PC3=as.numeric(topPC1[,3])
target$Seq.PC4=as.numeric(topPC1[,4])
target$Seq.PC5=as.numeric(topPC1[,5])
target$Seq.PC6=as.numeric(topPC1[,6])
target$Seq.PC7=as.numeric(topPC1[,7])
head(target)

# pre - pca/var ------------
PCAdatlist<-list()
PCAdatlist$seq$topvar<-seqtopvar
PCAdatlist$pre<-PCvar(normExpr)

plotvar(pclist = PCAdatlist, expPCs = "pre",
        outdir = figdir,fname = "screes.pdf")

# pre - heatmaps  --------------
plothms(pclist = PCAdatlist, expPCs = "pre",
        conds = conds, target = target,
        outdir = figdir, fname = "hms1.pdf")

# pairsdat ---------------
library(WGCNA)
cond=labels2colors(as.numeric(factor(target$GT))) 
n=dim(conds)[2]
pairsdat <- data.frame(
  GT = as.factor(target$GT),
  # RIN = as.numeric(target$RIN),
  Seq_PC1 = as.numeric(target$Seq.PC1),
  Seq_PC2 = as.numeric(target$Seq.PC2),
  Seq_PC3 = as.numeric(target$Seq.PC3),
  Seq_PC4 = as.numeric(target$Seq.PC4),
  Seq_PC5 = as.numeric(target$Seq.PC5),
  Seq_PC6 = as.numeric(target$Seq.PC6),
  Seq_PC7 = as.numeric(target$Seq.PC7)
)

# pre - cor matrix --------------------

cormatrix(pclist = PCAdatlist, expPCs = "pre", fname = "CorrPlot_Pre.pdf")


# correction setup -----------------
GT<-as.numeric(as.factor(target$GT))-1
seqPC1 <- as.numeric(target$Seq.PC1)
seqPC2 <- as.numeric(target$Seq.PC2)
seqPC3 <- as.numeric(target$Seq.PC3)
seqPC4 <- as.numeric(target$Seq.PC4)
seqPC5 <- as.numeric(target$Seq.PC5)
seqPC6 <- as.numeric(target$Seq.PC6)
seqPC7 <- as.numeric(target$Seq.PC7)



dat.cor<-list()
# cor1 spc1 -------------
modVars<-c("GT","seqPC1")
corname<-"SeqPC1"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr1.pdf")
mdsCompPlot(mds.list = mdslist,fname = "MDS_Comp1.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC1_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr1.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr1.pdf")



# cor2 spc2 -------------
modVars<-c("GT","seqPC2")
corname<-"SeqPC2"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr2.pdf")
mdsCompPlot(mds.list = mdslist,fname = "MDS_Comp2.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC2_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr2.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr2.pdf")


# cor2 spc3 -------------
modVars<-c("GT","seqPC3")
corname<-"SeqPC3"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr3.pdf")
mdsCompPlot(mds.list = mdslist,fname = "MDS_Comp3.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC3_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr3.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr3.pdf")

# cor4 spc12 -------------
modVars<-c("GT","seqPC1","seqPC2")
corname<-"SeqPC12"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr4.pdf")
mdsCompPlot(mds.list = mdslist[c("pre",corname)],fname = "MDS_Comp4.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC12_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr4.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr4.pdf")


# cor5 spc123 -------------
modVars<-c("GT","seqPC1","seqPC2","seqPC3")
corname<-"SeqPC123"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr5.pdf")
mdsCompPlot(mds.list = mdslist[c("pre",corname)],fname = "MDS_Comp5.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC123_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr5.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr5.pdf")



# cor6 spc1237 -------------
modVars<-c("GT","seqPC1","seqPC2","seqPC3","seqPC7")
corname<-"SeqPC1237"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr6.pdf")
mdsCompPlot(mds.list = mdslist[c("pre",corname)],fname = "MDS_Comp6.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC1237_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr6.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr6.pdf")




# cor7 spc1234 -------------
modVars<-c("GT","seqPC1","seqPC2","seqPC3","seqPC4")
corname<-"SeqPC1234"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr7.pdf")
mdsCompPlot(mds.list = mdslist[c("pre",corname)],fname = "MDS_Comp7.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC1234_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr7.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr7.pdf")


# cor8 spc12347 -------------
modVars<-c("GT","seqPC1","seqPC2","seqPC3","seqPC4","seqPC7")
corname<-"SeqPC12347"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr8.pdf")
mdsCompPlot(mds.list = mdslist[c("pre",corname)],fname = "MDS_Comp8.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC12347_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr8.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr8.pdf")


# dea------
datexprList<-dat.cor[c("SeqPC123","SeqPC1237","SeqPC1234","SeqPC12347")]
GT<-factor(target$GT)
design<-model.matrix(~0+GT)
colnames(design)<-levels(GT)
head(design)
contrasts<-makeContrasts(Hem - WT, levels=design)
contrasts
an<-fData(eset)
an$Gene.description<-gsub("[[:space:]]\\[.+","",an$Gene.description)
names(an)[names(an)=="Chromosome.scaffold.name"]<-"Chr"

fitList<-lapply(datexprList, function(x){
  fit <- lmFit(x, design)
  fit$genes<-an[rownames(x),c("Gene.name","Gene.description")]
  fit.cont<-contrasts.fit(fit, contrasts)
  fit.cont<-eBayes(fit.cont)
  return(fit.cont)
})



dtlist<-list()

for(cor.fit in names(fitList)){
  for(lfc in c(0,0.5,1,1.5)){
    for(pval in c(0.05,0.1)){
      dts<-as.data.frame(summary(decideTests(fitList[[cor.fit]],lfc = lfc,p.value = pval)))
      dts$lfc<-lfc
      dts$fdr<-pval
      dts$Correction<-cor.fit
      dtlist[[paste0(cor.fit,"lfc",lfc,"p",pval)]]<-dts
    }
  }
}
dtdat<-do.call(rbind,dtlist)
names(dtdat)[1:2]<-c("Dir","Cont")


dtdat$lfc<-factor(dtdat$lfc)
dtdat$fdr<-factor(dtdat$fdr)
dtdat$Dir<-factor(dtdat$Dir,levels = c("Up","NotSig","Down"))
dtdatnolfc<-dtdat[dtdat$lfc==0,]
ggplot(dtdatnolfc,aes(fdr,Freq,fill=Dir,label=Freq))+
  geom_col()+
  geom_label(position = position_stack(vjust = 0.5),show.legend = F)+
  facet_wrap(~Correction,ncol = 4)+
  theme_bw()+
  theme(strip.text = element_text(face = "bold"))

ggsave(paste0(figdir,"dea.pdf"),width = 8,height = 3)


# tt<-topTable(fit.cont,coef = "Hem - WT")
# tt2<-data.frame(Symb=tt$Gene.name,
#                 Name=tt$Gene.description,
#                 logFC=round(tt$logFC,2),
#                 P=tt$P.Value,
#                 Padj=tt$adj.P.Val)

ttlist<-lapply(fitList,function(x){
  tt<-topTable(x,coef = "Hem - WT")
  tt2<-data.frame(Symb=tt$Gene.name,
                  Name=tt$Gene.description,
                  logFC=round(tt$logFC,2),
                  P=tt$P.Value,
                  Padj=tt$adj.P.Val)
  return(tt2)
})

library(knitr)


for(cor.fit in names(ttlist)){
  sink(file = paste0(outname,"/tex/tt",cor.fit,".tex"))
  cat("\\begin{table}\n")
  cat("\\begin{adjustbox}{width=0.9\\textwidth, totalheight=\\textheight-2\\baselineskip-2\\baselineskip,keepaspectratio}\n")
  print(kable(ttlist[[cor.fit]],format = "latex",row.names = F))
  cat("\\end{adjustbox}\n")
  cat("\\caption{",cor.fit,"}",sep = "")
  cat("\n")
  cat("\\end{table}\n")
  sink()
}



ttlistfull<-lapply(fitList,function(x){
  tt<-topTable(x,coef = "Hem - WT",adjust.method = "BH",number = "inf")
  return(tt)
})

saveRDS(ttlistfull, file = paste0(datdir,"ttlist.rds"))

tt<-topTable(fitList$SeqPC12347,coef = "Hem - WT",adjust.method = "BH",number = "inf")
saveRDS(tt, file = paste0(datdir,"tt_seqPC12347.rds"))
write.csv(tt, file = paste0(datdir,"tt_seqPC12347.csv"))

eset.cor<-filtNormEset
dim(pData(eset.cor))
dim(eset.cor)
dim(dat.cor$SeqPC12347)
exprs(eset.cor)<-as.matrix(dat.cor$SeqPC12347)

eset<-readRDS("data/rawCountEsetiba.rds")
dim(eset)
ibadat<-pData(eset[,rownames(pData(eset.cor))])
dim(ibadat)
match(rownames(pData(eset.cor)),rownames(ibadat))
pData(eset.cor)<-ibadat
an<-fData(eset.cor)
an$Gene.description<-gsub("[[:space:]]\\[.+","",an$Gene.description)
fData(eset.cor)<-an

saveRDS(eset.cor,"dc1m_spc12347eset.rds")


save.image(file = paste0(datdir,"rdata.rda"))


# load(paste0(datdir,"rdata.rda"))
