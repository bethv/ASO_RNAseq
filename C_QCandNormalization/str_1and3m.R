rm(list = ls())

source("QCfunctions.R")
grp="str"
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
outliers=c("7648_dc","7675_dc","7609_dc","7664_dc")

eset<-eset[,!rownames(pData(eset))%in%outliers]
eset<-eset[,eset$tissue==grp]
dim(eset);table(eset$GTtime)

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
                              design = ~ GTtime)
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
# filtNormEset<-readRDS(filtNormEsetFile)

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
conds<-conds[,c("sample","RIN","GT","Time","GTtime")]
conds2<-conds[,c("sample","GTtime")]




mdslist<-list()
mdslist$pre<-makemds.set(data = normDat)

mdsplot1v2(mdslist$pre,figdir,"MDS_pre.pdf")
mdsplotAll(mdslist$pre,figdir,"MDS_pre_all.pdf")


# rin plot -------------
# library(tidyverse)
# library(gridExtra)
# library(grid)

# rinplot<-ggplot(conds, aes(Time,RIN,fill=GT))+
#   geom_boxplot(position = position_dodge(1))+
#   geom_dotplot(binaxis = "y",stackdir = "center",position = position_dodge(1))+
#   theme_bw()+
#   ggtitle(grp)

# tts<-list()

# tts$mainGT<-t.test(conds$RIN[conds$GT=="Hem"],
#           conds$RIN[conds$GT=="WT"])

# tts$mainTime<-t.test(conds$RIN[conds$Time==1],
#                    conds$RIN[conds$Time==3])

# tts$GT1<-t.test(conds$RIN[conds$GT=="Hem"&conds$Time==1],
#                    conds$RIN[conds$GT=="WT"&conds$Time==1])

# tts$GT3<-t.test(conds$RIN[conds$GT=="Hem"&conds$Time==3],
#                 conds$RIN[conds$GT=="WT"&conds$Time==3])

# pvals<-lapply(tts, function(x){
#   p<-round(x$p.value,3)
#   return(p)
# })


# ptab<-as.data.frame(cbind(c("Main effect of GT",
#               "Main effect of time",
#               "GT at 1 month",
#               "GT at 3 months"),
#             do.call(rbind,pvals)))

# names(ptab)<-c("Comparison","p val")

# rinfig<-arrangeGrob(rinplot, tableGrob(ptab,rows = NULL), ncol=2)

# ggsave(plot = rinfig,paste0(figdir,"rintt.pdf"),width = 7.5,height = 4)

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
cond=labels2colors(as.numeric(factor(target$GTtime))) 
n=dim(conds)[2]
t1m<-target$Time==1
GT_1m<-target$GT
GT_1m[!t1m]<-NA
GT_3m<-target$GT
GT_3m[t1m]<-NA

pairsdat <- data.frame(
  GT = as.factor(target$GT),
  Time = as.factor(target$Time),
#   GTtime = as.factor(target$GTtime),
  GT_1m = as.factor(GT_1m),
  GT_3m = as.factor(GT_3m),  
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


# # pre explore interaction-----------
# GT<-as.numeric(as.factor(target$GT))-1
# Time<-as.numeric(as.factor(target$Time))-1

# exp.pcdat<-cbind.data.frame(GT=as.numeric(as.factor(target$GT))-1,
#                             Time=as.numeric(as.factor(target$Time))-1,
#                             PCAdatlist$pre$TopPC1)
# names(exp.pcdat)[3:7]<-paste0("ExpPC",seq(1,5,1))

# fit<-lm(ExpPC5~GT+Time+GT*Time, data=exp.pcdat);summary(fit)


# correction setup -----------------
GT<-as.numeric(as.factor(target$GT))-1
Time<-as.numeric(as.factor(target$Time))-1
# RIN <- as.numeric(target$RIN)
seqPC1 <- as.numeric(target$Seq.PC1)
seqPC2 <- as.numeric(target$Seq.PC2)
seqPC3 <- as.numeric(target$Seq.PC3)
seqPC4 <- as.numeric(target$Seq.PC4)
seqPC5 <- as.numeric(target$Seq.PC5)
seqPC6 <- as.numeric(target$Seq.PC6)
seqPC7 <- as.numeric(target$Seq.PC7)

dat.cor<-list()


# cor1 spc1 -------------
modVars<-c("GT","Time","seqPC1")
corname<-"SeqPC1"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr1.pdf")
mdsCompPlot(mds.list = mdslist[c("pre",corname)],fname = "MDS_Comp1.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC1_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr1.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr1.pdf")

# cor2 spc12 -------------
modVars<-c("GT","Time","seqPC1","seqPC2")
corname<-"SeqPC12"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr2.pdf")
mdsCompPlot(mds.list = mdslist[c("pre","SeqPC1",corname)],fname = "MDS_Comp2.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC12_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr2.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr2.pdf")

# cor3 spc123 -------------
modVars<-c("GT","Time","seqPC1","seqPC2","seqPC3")
corname<-"SeqPC123"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr3.pdf")
mdsCompPlot(mds.list = mdslist[c("pre","SeqPC12",corname)],fname = "MDS_Comp3.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC123_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr3.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr3.pdf")



# cor4 spc1234 -------------
modVars<-c("GT","Time","seqPC1","seqPC2","seqPC3","seqPC4")
corname<-"SeqPC1234"
dat.cor[[corname]]<-correct()

mdslist[[corname]]<-makemds.set(data = dat.cor[[corname]])
mdsplot1v2(mds.set = mdslist[[corname]],fname = "MDS_Corr4.pdf")
mdsCompPlot(mds.list = mdslist[c("pre","SeqPC123",corname)],fname = "MDS_Comp4.pdf")
mdsplotAll(mds.set = mdslist[[corname]],fname = "MDS_SeqPC1234_all.pdf")

PCAdatlist[[corname]]<-PCvar(dat.cor[[corname]])
plothms(pclist = PCAdatlist,
        expPCs = corname,
        fname = "hmsCorr4.pdf")
cormatrix(pclist = PCAdatlist,
          expPCs = corname,
          fname = "CorrPlot_Corr4.pdf")


# dea------

datexprList<-dat.cor[c("SeqPC1","SeqPC12","SeqPC123","SeqPC1234")]
GTtime<-factor(target$GTtime)
design<-model.matrix(~0+GTtime)
colnames(design)<-levels(GTtime)
head(design)
an<-fData(eset)
an$Gene.description<-gsub("[[:space:]]\\[.+","",an$Gene.description)
names(an)[names(an)=="Chromosome.scaffold.name"]<-"Chr"

contrasts<-makeContrasts(Hem_1 - WT_1,Hem_3 - WT_3, levels=design)
contrasts

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

dtdat$Correction<-factor(dtdat$Correction, levels=c("SeqPC1","SeqPC12","SeqPC123","SeqPC1234"))

dtdat$lfc<-factor(dtdat$lfc)
dtdat$fdr<-factor(dtdat$fdr)
dtdat$Dir<-factor(dtdat$Dir,levels = c("Up","NotSig","Down"))
dtdatnolfc<-dtdat[dtdat$lfc==0,]
ggplot(dtdatnolfc,aes(fdr,Freq,fill=Dir,label=Freq))+
  geom_col()+
  geom_label(position = position_stack(vjust = 0.5),show.legend = F)+
  facet_grid(Cont~Correction)+
  theme_bw()+
  theme(strip.text = element_text(face = "bold"))

ggsave(paste0(figdir,"dea.pdf"),width = 8,height = 5)


ttlist1m<-lapply(fitList,function(x){
  tt<-topTable(x,coef = "Hem_1 - WT_1")
  tt2<-data.frame(Symb=tt$Gene.name,
                  Name=tt$Gene.description,
                  logFC=round(tt$logFC,2),
                  P=tt$P.Value,
                  Padj=tt$adj.P.Val)
  return(tt2)
})


ttlist3m<-lapply(fitList,function(x){
  tt<-topTable(x,coef = "Hem_3 - WT_3")
  tt2<-data.frame(Symb=tt$Gene.name,
                  Name=tt$Gene.description,
                  logFC=round(tt$logFC,2),
                  P=tt$P.Value,
                  Padj=tt$adj.P.Val)
  return(tt2)
})


library(knitr)


for(cor.fit in names(ttlist1m)){
  sink(file = paste0(outname,"/tex/tt1m",cor.fit,".tex"))
  cat("\\begin{table}\n")
  cat("\\begin{adjustbox}{width=0.9\\textwidth, totalheight=\\textheight-2\\baselineskip-2\\baselineskip,keepaspectratio}\n")
  print(kable(ttlist1m[[cor.fit]],format = "latex",row.names = F))
  cat("\\end{adjustbox}\n")
  cat("\\caption{",cor.fit,"}",sep = "")
  cat("\n")
  cat("\\end{table}\n")
  sink()
}


for(cor.fit in names(ttlist3m)){
  sink(file = paste0(outname,"/tex/tt3m",cor.fit,".tex"))
  cat("\\begin{table}\n")
  cat("\\begin{adjustbox}{width=0.9\\textwidth, totalheight=\\textheight-2\\baselineskip-2\\baselineskip,keepaspectratio}\n")
  print(kable(ttlist3m[[cor.fit]],format = "latex",row.names = F))
  cat("\\end{adjustbox}\n")
  cat("\\caption{",cor.fit,"}",sep = "")
  cat("\n")
  cat("\\end{table}\n")
  sink()
}



save.image(file = paste0(datdir,"rdata.rda"))



# # abort and try 1 and 3 m sep-------


# rm(list = ls())
# # time sep - setup ---------
# source("QCfunctions.R")
# grp="str"
# datdir<-paste0(grp,"/data/")
# figdir<-paste0(grp,"/figs/")
# eset<-readRDS("data/rawCountEset.rds") 
# eset$Reads<-as.numeric(gsub(",","",eset$Reads))
# outliers=c("7648_dc","7675_dc")
# eset<-eset[,!rownames(pData(eset))%in%outliers]
# eset<-eset[,eset$tissue==grp]
# times<-c(1,3)
# tdatdirs<-paste0(datdir,grp,times,"m/")
# tfigdirs<-paste0(figdir,grp,times,"m/")
# lapply(tdatdirs,function(x) dir.create(x,showWarnings = FALSE))
# lapply(tfigdirs,function(x) dir.create(x,showWarnings = FALSE))
# tgrps<-paste0(grp,times,"m")
# esetList<-lapply(times,function(x){
#   return(eset[,eset$Time==x])
# })

# lapply(esetList,dim)

# # byT - filter/vsd ---------------

# filteredEsetList<-list()
# cdsList<-list()
# normDatList<-list()
# filtNormEsetList<-list()

# for(a in 1:2){
#   datdir<-tdatdirs[[a]]
  
#   eset<-esetList[[a]]
#   temp<-as.data.frame(exprs(eset))
#   temp$flag<-0
#   for (i in 1:ncol(exprs(eset))){
#     temp$flag<-temp$flag+(temp[,i]>=10)
#   }
#   fiftyPercent=round(.5*dim(eset)[2],digits = 0)
#   fiftyPercent
#   table(temp$flag>=fiftyPercent) 
#   keep=temp$flag>=fiftyPercent
#   filteredEset<-eset[keep,]
#   dim(filteredEset)
#   filtEsetFile<-paste0(datdir,"filteredEset.rds")
#   saveRDS(filteredEset,filtEsetFile)
#   filteredEsetList[[a]]<-filteredEset
#   cds <- DESeqDataSetFromMatrix(countData = exprs(filteredEset),
#                                 colData = pData(filteredEset),
#                                 design = ~ GT)
#   cds <- estimateSizeFactors(cds)
#   cds <- estimateDispersions(cds)
#   cdsList[[a]]<-cds
#   vsd <- getVarianceStabilizedData(cds)
#   normDat<-vsd
#   saveRDS(normDat,paste0(datdir,"NormCounts.rds"))
#   normDatList[[a]]<-normDat
#   #make into eset
#   filtNormEset<-filteredEset
#   exprs(filtNormEset)<-normDat
#   filtNormEsetFile<-paste0(datdir,"filtNormEset.rds")
#   saveRDS(filtNormEset,filtNormEsetFile)
#   filtNormEsetList[[a]]<-filtNormEset
# }


# # byT - density plots -------------------

# for(a in 1:2){
#   figdir<-tfigdirs[[a]]
#   eset<-esetList[[a]]
#   filteredEset<-filteredEsetList[[a]]
#   filtNormEset<-filtNormEsetList[[a]]
#   cds<-cdsList[[a]]
#   pdf(paste0(figdir,"densityplots.pdf"),width = 7.5*2, height = 9*2)
  
#   lm=rbind(c(13,8,8,8),
#            c(9,1,2,2),
#            c(10,3,4,4),
#            c(11,5,6,6),
#            c(12,7,7,7))
  
#   layout(mat = lm,
#          heights = c(1,6,6,6,8))
  
#   dlines(exprs(eset), main = paste0("n=",dim(eset)[1]))
#   dboxes(exprs(eset))
#   dlines(exprs(filteredEset),main = paste0("n=",dim(filteredEset)[1]))
#   dboxes(exprs(filteredEset))
#   dlines(exprs(filtNormEset),log = F,main = paste0("n=",dim(filtNormEset)[1]))
#   dboxes(exprs(filtNormEset),log = F)
#   plotDispEsts(cds)
  
#   par(mar = c(0,0,0,0))
#   plot.new()
#   text(0.5,0.5,paste0("striatum ",times[a],"m, n=",dim(eset)[2]),cex=3)
  
#   plot.new();text(0.5,0.5,"raw",cex=3)
#   plot.new();text(0.5,0.5,"filtered",cex=3)
#   plot.new();text(0.5,0.5,"VSD",cex=3)
#   dev.off()
  
#   pdf(paste0(figdir,"densityplots2.pdf"),width = 8, height = 9)
#   lm=rbind(c(7,1,2,2),
#            c(8,3,4,4),
#            c(9,5,6,6))
#   layout(lm)
#   dlines(exprs(eset), main = paste0("n=",dim(eset)[1]))
#   dboxes(exprs(eset))
#   dlines(exprs(filteredEset),main = paste0("n=",dim(filteredEset)[1]))
#   dboxes(exprs(filteredEset))
#   dlines(exprs(filtNormEset),log = F,main = paste0("n=",dim(filtNormEset)[1]))
#   dboxes(exprs(filtNormEset),log = F)
  
#   par(mar = c(0,0,0,0))
#   plot.new();text(0.5,0.5,"raw",cex=3)
#   plot.new();text(0.5,0.5,"filtered",cex=3)
#   plot.new();text(0.5,0.5,"VSD",cex=3)
#   dev.off()
  
#   pdf(paste0(figdir,"meanvar.pdf"),width = 5, height = 5)
#   plotDispEsts(cds)
#   dev.off()
  
# }



# # byT - pre - mds 1v2 ----------

# condsList<-list()
# conds2List<-list()
# mdsLists<-list(str1m=list(),str3m=list())
# for(a in 1:2){
#   figdir<-tfigdirs[[a]]
#   datdir<-tdatdirs[[a]]
#   filtNormEset<-filtNormEsetList[[a]]
#   mdslist<-mdsLists[[a]]
#   normDat<-normDatList[[a]]
#   conds<-pData(filtNormEset)
#   conds$sample<-rownames(conds)
#   conds$RIN<-as.numeric(conds$RIN)
#   conds<-conds[,c("sample","RIN","GT")]
#   conds2<-conds[,c("sample","GT")]
#   condsList[[a]]<-conds2
#   conds2List[[a]]<-conds2
  
#   mdslist$pre<-list(
#     All=getplotdata(data = normDat, conds2 = conds2, kfit = 6,topn = NULL),
#     top500=getplotdata(data = normDat, conds2 = conds2, kfit = 6,topn = 500),
#     top1000=getplotdata(data = normDat, conds2 = conds2, kfit = 6,topn = 1000)
#   )
  
  
#   pdf(paste0(figdir,"MDS_pre.pdf"),width = 10, height = 5)
#   sizeText=1; par(mfrow=c(1,3))
#   for(i in 1:length(mdslist$pre)){
#     plotdat<-mdslist$pre[[i]]
#     makeplot(plotdat = plotdat,1,2)
#     # makeplot(plotdat = plotdat,3,4)
#     # makeplot(plotdat = plotdat,5,6)
#   }
#   dev.off()
  
#   # pre - mds all ----------
#   pdf(paste0(figdir,"MDS_pre_all.pdf"),width = 10, height = 5)
#   sizeText=.75; par(mar = c(0,0,0,0))
#   layout(mat = matrix(seq(1,16),ncol = 4, byrow = T),
#          heights = c(.5,4,4,4),widths = c(2,4,4,4))
#   plot.new()
#   plot.new();text(0.5,0.5,"All",cex=2)
#   plot.new();text(0.5,0.5,"Top 500 (var)",cex=2)
#   plot.new();text(0.5,0.5,"Top 1000 (var)",cex=2)
#   plot.new();text(0.5,0.5,"1 v 2",cex=2)
#   for(i in 1:length(mdslist$pre)){
#     plotdat<-mdslist$pre[[i]]
#     makeplot(plotdat = plotdat,1,2,main="",axes = F)
#   }
#   plot.new();text(0.5,0.5,"3 v 4",cex=2)
#   for(i in 1:length(mdslist$pre)){
#     plotdat<-mdslist$pre[[i]]
#     makeplot(plotdat = plotdat,3,4,main="",axes = F)
#   }
#   plot.new();text(0.5,0.5,"5 v 6",cex=2)
#   for(i in 1:length(mdslist$pre)){
#     plotdat<-mdslist$pre[[i]]
#     makeplot(plotdat = plotdat,5,6,main="",axes = F)
#   }
#   dev.off()
#   mdsLists[[a]]<-mdslist
# }






# # byT - picard --------------
# picard<-readRDS("data/picardTable.rds")
# picard<-picard[,c(
#   "N_unmapped", "N_multimapping","N_noFeature","N_ambiguous", #star count info
  
#   #AlignmentSummaryMetrics
#   "PF_READS", "PF_HQ_ALIGNED_READS",
#   "PF_INDEL_RATE",  "PF_HQ_ERROR_RATE", 
#   "PF_ALIGNED_BASES", "PF_HQ_ALIGNED_BASES", "PF_HQ_ALIGNED_Q20_BASES", 
#   "PF_MISMATCH_RATE", 
  
#   #rnaseqmetrics
#   "PF_BASES", "CODING_BASES", "PCT_CODING_BASES", 
#   "UTR_BASES", "PCT_UTR_BASES",
#   "INTRONIC_BASES", "PCT_INTRONIC_BASES",
#   "INTERGENIC_BASES", "PCT_INTERGENIC_BASES", 
#   "CORRECT_STRAND_READS", "PCT_CORRECT_STRAND_READS",
#   "INCORRECT_STRAND_READS",
#   "PCT_MRNA_BASES", "PCT_USABLE_BASES",
#   "MEDIAN_CV_COVERAGE", 
  
#   #DuplicationMetrics
#   "UNPAIRED_READ_DUPLICATES",  "PERCENT_DUPLICATION")]

# targetList<-list()
# seqtopvarList<-list()
# normExprList<-list()
# for(a in 1:2){
#   conds<-condsList[[a]]
#   filtNormEset<-filtNormEsetList[[a]]
#   picardgrp<-picard[rownames(conds),]
#   picardgrp<-as.data.frame(apply(picardgrp, 2, as.numeric))
#   rownames(picardgrp)<-rownames(conds);head(picardgrp)
#   target<-cbind(conds,picardgrp);dim(target)
#   n=dim(conds)[2]
#   target<-target[,!(colSums(is.na(target)) > 0)];dim(target)
#   temp1<-target[,-c(1:n)];print(dim(temp1))
#   temp1<-temp1[,colVars(as.matrix(temp1))!=0];dim(temp1)
#   temp1<-t(scale(temp1))
#   temp1<-as.data.frame(t(temp1))
#   target<-cbind(target[,1:n],temp1)

#   normExpr<-as.data.frame(exprs(filtNormEset))
#   normExprList[[a]]<-normExpr
#   unique(match(names(normExpr),target$sample)==seq(1,nrow(target),1))
#   tar1=target[,c((1+dim(conds)[2]):dim(target)[2])]
#   thisdat <- t(scale(tar1,scale=F))
#   PC.metadata <- prcomp(thisdat,center=F);
#   topPC1 <- PC.metadata$rotation[,1:10];
#   varexp <- (PC.metadata$sdev)^2 / sum(PC.metadata$sdev^2)
#   seqtopvar <- varexp[1:10]
#   colnames(topPC1) <- paste("Seq.PC\n",colnames(topPC1)," (",signif(100*seqtopvar[1:7],2),"%)",sep="")
#   target$Seq.PC1=as.numeric(topPC1[,1])
#   target$Seq.PC2=as.numeric(topPC1[,2])
#   target$Seq.PC3=as.numeric(topPC1[,3])
#   target$Seq.PC4=as.numeric(topPC1[,4])
#   target$Seq.PC5=as.numeric(topPC1[,5])
#   target$Seq.PC6=as.numeric(topPC1[,6])
#   target$Seq.PC7=as.numeric(topPC1[,7])
#   targetList[[a]]<-target
#   seqtopvarList[[a]]<-seqtopvar
  
# }


# # pre - pca/var,heatmaps ------------
# PCAdatLists<-list(str1m=list(), str3m=list())
# r2matLists<-list(str1m=list(), str3m=list())

# for(a in 1:2){
#   normExpr<-normExprList[[a]]
#   figdir<-tfigdirs[[a]]
#   conds<-condsList[[a]]
#   target<-targetList[[a]]
#   PCAdatLists[[a]][["pre"]]<-PCvar(normExpr)
#   pdf(paste0(figdir,"screes.pdf"),height = 4,width = 7.5)
#   par(mfrow=c(1,2))
#   plot(PCAdatLists[[a]][["pre"]][["topvar"]], xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", 
#        ylab = "variance explained", main="Variance explained by Expression PCs",cex.main=1)
#   plot(seqtopvarList[[a]], xlim = c(0, 10), type = "b", pch = 16, xlab = "principal components", 
#        ylab = "variance explained", main="Variance explained by sequencing PCs",cex.main=1)
#   dev.off()
  
#   r2matLists[[a]][["pre"]]<-get.r2mats(
#     PCAdatLists[[a]][["pre"]][["TopPC1"]],
#     conds = conds,target = target)
#   hm1<-pheatmap(r2matLists[[a]][["pre"]][["r2mat_meta"]],
#                 cluster_cols = F,cluster_rows=FALSE, fontsize = 8, 
#                 silent=T,display_numbers = TRUE,fontsize_number = 8,
#                 cellwidth = 20, cellheight = 20)
#   hm2<-pheatmap(t(r2matLists[[a]][["pre"]][["r2mat_seq"]]),
#                 clustering_method="average",fontsize = 8,silent = T,
#                 cluster_cols=FALSE, display_numbers = TRUE,
#                 fontsize_number = 8, cellwidth = 40)
#   hms<-arrangeGrob(grobs = list(hm1[[4]],hm2[[4]]), 
#                    nrow = 1,widths = c(3,6))
#   # hms<-grid.arrange(hms)
#   ggsave(plot = hms,filename = paste0(figdir,"hms1.pdf"),
#          width = 8,height = 5)
#   }




# # byT - pairsdat ---------------
# pairsdatList<-list()
# condList<-list()
# for(a in 1:2){
#   target<-targetList[[a]]
#   cond=labels2colors(as.numeric(factor(target$GT))) 
#   n=dim(conds)[2]
#   pairsdat <- data.frame(
#     GT = as.factor(target$GT),
#     # RIN = as.numeric(target$RIN),
#     Seq_PC1 = as.numeric(target$Seq.PC1),
#     Seq_PC2 = as.numeric(target$Seq.PC2),
#     Seq_PC3 = as.numeric(target$Seq.PC3),
#     Seq_PC4 = as.numeric(target$Seq.PC4),
#     Seq_PC5 = as.numeric(target$Seq.PC5),
#     Seq_PC6 = as.numeric(target$Seq.PC6),
#     Seq_PC7 = as.numeric(target$Seq.PC7)
#   )
#   pairsdatList[[a]]<-pairsdat
#   condList[[a]]<-cond
# }


# # byT - pre - cor matrix --------------------

# for(a in 1:2){
#   pdf(paste0(tfigdirs[[a]],"CorrPlot_PreN.pdf"),width = 9, height=8)
#   pairs(cbind(PCAdatLists[[a]][["pre"]][["TopPC1"]],
#               pairsdatList[[a]]),
#         col= condList[[a]],pch=19,upper.panel = panel.cor, cex.labels=1,gap = .25)
#   dev.off()
  
# }






# # 1m - cor1 spc4 -------------
# target<-targetList[[1]]
# normExpr<-normExprList[[1]]
# GT<-as.numeric(as.factor(target$GT))-1
# seqPC1 <- as.numeric(target$Seq.PC1)
# seqPC2 <- as.numeric(target$Seq.PC2)
# seqPC3 <- as.numeric(target$Seq.PC3)
# seqPC4 <- as.numeric(target$Seq.PC4)
# seqPC5 <- as.numeric(target$Seq.PC5)
# seqPC6 <- as.numeric(target$Seq.PC6)
# seqPC7 <- as.numeric(target$Seq.PC7)

# X =model.matrix(~GT+seqPC4)
# head(X)
# regColIdx=3
# Y = normExpr
# beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
# b = as.data.frame(t(beta))
# to_regress = (as.matrix(X[,regColIdx:dim(X)[2],drop=FALSE]) %*% 
#                 (as.matrix(beta[regColIdx:dim(X)[2],,drop=FALSE])))
# datExprfit = normExpr - t(to_regress) 

# dat.cor1m<-list()
# dat.cor1m$SeqPC4<-datExprfit

# # 1m cor1 - mds 1v2 ---------------
# mdslist1m<-mdsLists[[1]]
# conds2<-conds2List[[1]]
# mdslist1m$SeqPC4<-list(
#   All=getplotdata(data = dat.cor1m$SeqPC4, conds2 = conds2, kfit = 6,topn = NULL),
#   top500=getplotdata(data = dat.cor1m$SeqPC4, conds2 = conds2, kfit = 6,topn = 500),
#   top1000=getplotdata(data = dat.cor1m$SeqPC4, conds2 = conds2, kfit = 6,topn = 1000)
# )

# figdir<-tfigdirs[[1]]
# pdf(paste0(figdir,"MDS_Corr1.pdf"),width = 10, height = 5)
# sizeText=1; par(mfrow=c(1,3))
# for(i in 1:length(mdslist1m$SeqPC4)){
#   plotdat<-mdslist1m$SeqPC4[[i]]
#   makeplot(plotdat = plotdat,1,2)
# }
# dev.off()

# # 1m cor1 - mds comp ---------------
# pdf(paste0(figdir,"MDS_Comp1.pdf"),width = 10, height = 5)
# sizeText=.75; par(mar = c(0,0,0,0))
# layout(mat = matrix(seq(1,12),ncol = 4, byrow = T),
#        heights = c(.5,4,4),widths = c(3,4,4,4))
# plot.new();text(0.5,0.5,"Coord 1 v 2",cex=2)
# plot.new();text(0.5,0.5,"All",cex=2)
# plot.new();text(0.5,0.5,"Top 500 (var)",cex=2)
# plot.new();text(0.5,0.5,"Top 1000 (var)",cex=2)
# plot.new();text(0.5,0.5,"Pre",cex=2)
# for(i in 1:length(mdslist1m$pre)){
#   plotdat<-mdslist1m$pre[[i]]
#   makeplot(plotdat = plotdat,1,2,main="",axes = F)
# }
# plot.new();text(0.5,0.5,"SeqPC4",cex=2)
# for(i in 1:length(mdslist1m$SeqPC4)){
#   plotdat<-mdslist1m$SeqPC4[[i]]
#   makeplot(plotdat = plotdat,1,2,main="",axes = F)
# }
# dev.off()

# # 1m cor1 - mds all ---------------
# pdf(paste0(figdir,"MDS_SeqPC1_all.pdf"),width = 10, height = 5)
# sizeText=.75; par(mar = c(0,0,0,0))
# layout(mat = matrix(seq(1,16),ncol = 4, byrow = T),
#        heights = c(.5,4,4,4),widths = c(2,4,4,4))
# plot.new()
# plot.new();text(0.5,0.5,"All",cex=2)
# plot.new();text(0.5,0.5,"Top 500 (var)",cex=2)
# plot.new();text(0.5,0.5,"Top 1000 (var)",cex=2)
# plot.new();text(0.5,0.5,"1 v 2",cex=2)
# for(i in 1:length(mdslist1m$SeqPC4)){
#   plotdat<-mdslist1m$SeqPC4[[i]]
#   makeplot(plotdat = plotdat,1,2,main="",axes = F)
# }
# plot.new();text(0.5,0.5,"3 v 4",cex=2)
# for(i in 1:length(mdslist1m$SeqPC4)){
#   plotdat<-mdslist1m$SeqPC4[[i]]
#   makeplot(plotdat = plotdat,3,4,main="",axes = F)
# }
# plot.new();text(0.5,0.5,"5 v 6",cex=2)
# for(i in 1:length(mdslist1m$SeqPC4)){
#   plotdat<-mdslist1m$SeqPC4[[i]]
#   makeplot(plotdat = plotdat,5,6,main="",axes = F)
# }
# dev.off()

# # 1m cor1 - hms ---------------
# PCAdatlist1m<-PCAdatLists[[1]]
# r2matlist1m<-r2matLists[[1]]
# PCAdatlist1m$SeqPC4<-PCvar(dat.cor1m$SeqPC4)
# r2matlist1m$SeqPC4<-get.r2mats(topPC = PCAdatlist1m$SeqPC4$TopPC1,conds = condsList[[1]], target = targetList[[1]])

# hm1<-pheatmap(r2matlist1m$SeqPC4$r2mat_meta,cluster_cols = F,cluster_rows=FALSE, fontsize = 8, 
#               silent=T,display_numbers = TRUE,fontsize_number = 8,cellwidth = 20, cellheight = 20)
# hm2<-pheatmap(t(r2matlist1m$SeqPC4$r2mat_seq),clustering_method="average",fontsize = 8,silent = T,
#               cluster_cols=FALSE, display_numbers = TRUE,fontsize_number = 8, cellwidth = 40)
# hms<-arrangeGrob(grobs = list(hm1[[4]],hm2[[4]]), 
#                  nrow = 1,widths = c(3,6))
# ggsave(plot = hms,filename = paste0(figdir,"hmsCorr1.pdf"),
#        width = 8,height = 5)


# # 1m cor1 - matrix --------------------
# pdf(paste0(figdir,"CorrPlot_Corr1N.pdf"),width = 9, height=8)
# pairs(cbind(PCAdatlist1m$SeqPC4$TopPC1,pairsdatList[[1]]),
#       col= condList[[1]],pch=19,upper.panel = panel.cor, cex.labels=1,gap = .25)
# dev.off()

















# # 3m - cor1 spc2 -------------
# target<-targetList[[2]]
# normExpr<-normExprList[[2]]
# GT<-as.numeric(as.factor(target$GT))-1
# seqPC1 <- as.numeric(target$Seq.PC1)
# seqPC2 <- as.numeric(target$Seq.PC2)
# seqPC3 <- as.numeric(target$Seq.PC3)
# seqPC4 <- as.numeric(target$Seq.PC4)
# seqPC5 <- as.numeric(target$Seq.PC5)
# seqPC6 <- as.numeric(target$Seq.PC6)
# seqPC7 <- as.numeric(target$Seq.PC7)

# X =model.matrix(~GT+seqPC2)
# head(X)
# regColIdx=3
# Y = normExpr
# beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
# b = as.data.frame(t(beta))
# to_regress = (as.matrix(X[,regColIdx:dim(X)[2],drop=FALSE]) %*% 
#                 (as.matrix(beta[regColIdx:dim(X)[2],,drop=FALSE])))
# datExprfit = normExpr - t(to_regress) 

# dat.cor3m<-list()
# dat.cor3m$SeqPC2<-datExprfit

# # 3m cor1 - mds 1v2 ---------------
# mdslist3m<-mdsLists[[2]]
# conds2<-conds2List[[2]]
# mdslist3m$SeqPC2<-list(
#     All=getplotdata(data = dat.cor3m$SeqPC2, conds2 = conds2, kfit = 6,topn = NULL),
#   top500=getplotdata(data = dat.cor3m$SeqPC2, conds2 = conds2, kfit = 6,topn = 500),
#   top1000=getplotdata(data = dat.cor3m$SeqPC2, conds2 = conds2, kfit = 6,topn = 1000)
# )

# figdir<-tfigdirs[[2]]
# pdf(paste0(figdir,"MDS_Corr1.pdf"),width = 10, height = 5)
# sizeText=1; par(mfrow=c(1,3))
# for(i in 1:length(mdslist3m$SeqPC2)){
#   plotdat<-mdslist3m$SeqPC2[[i]]
#   makeplot(plotdat = plotdat,1,2)
# }
# dev.off()

# # 3m cor1 - mds comp ---------------
# pdf(paste0(figdir,"MDS_Comp1.pdf"),width = 10, height = 5)
# sizeText=.75; par(mar = c(0,0,0,0))
# layout(mat = matrix(seq(1,12),ncol = 4, byrow = T),
#        heights = c(.5,4,4),widths = c(3,4,4,4))
# plot.new();text(0.5,0.5,"Coord 1 v 2",cex=2)
# plot.new();text(0.5,0.5,"All",cex=2)
# plot.new();text(0.5,0.5,"Top 500 (var)",cex=2)
# plot.new();text(0.5,0.5,"Top 1000 (var)",cex=2)
# plot.new();text(0.5,0.5,"Pre",cex=2)
# for(i in 1:length(mdslist3m$pre)){
#   plotdat<-mdslist3m$pre[[i]]
#   makeplot(plotdat = plotdat,1,2,main="",axes = F)
# }
# plot.new();text(0.5,0.5,"SeqPC2",cex=2)
# for(i in 1:length(mdslist3m$SeqPC2)){
#   plotdat<-mdslist3m$SeqPC2[[i]]
#   makeplot(plotdat = plotdat,1,2,main="",axes = F)
# }
# dev.off()

# # 3m cor1 - mds all ---------------
# pdf(paste0(figdir,"MDS_SeqPC2_all.pdf"),width = 10, height = 5)
# sizeText=.75; par(mar = c(0,0,0,0))
# layout(mat = matrix(seq(1,16),ncol = 4, byrow = T),
#        heights = c(.5,4,4,4),widths = c(2,4,4,4))
# plot.new()
# plot.new();text(0.5,0.5,"All",cex=2)
# plot.new();text(0.5,0.5,"Top 500 (var)",cex=2)
# plot.new();text(0.5,0.5,"Top 1000 (var)",cex=2)
# plot.new();text(0.5,0.5,"1 v 2",cex=2)
# for(i in 1:length(mdslist3m$SeqPC2)){
#   plotdat<-mdslist3m$SeqPC2[[i]]
#   makeplot(plotdat = plotdat,1,2,main="",axes = F)
# }
# plot.new();text(0.5,0.5,"3 v 4",cex=2)
# for(i in 1:length(mdslist1m$SeqPC4)){
#   plotdat<-mdslist3m$SeqPC2[[i]]
#   makeplot(plotdat = plotdat,3,4,main="",axes = F)
# }
# plot.new();text(0.5,0.5,"5 v 6",cex=2)
# for(i in 1:length(mdslist1m$SeqPC4)){
#   plotdat<-mdslist3m$SeqPC2[[i]]
#   makeplot(plotdat = plotdat,5,6,main="",axes = F)
# }
# dev.off()

# # 3m cor1 - hms ---------------
# PCAdatlist3m<-PCAdatLists[[2]]
# r2matlist3m<-r2matLists[[2]]
# PCAdatlist3m$SeqPC2<-PCvar(dat.cor3m$SeqPC2)
# r2matlist3m$SeqPC2<-get.r2mats(topPC = PCAdatlist3m$SeqPC2$TopPC1,
#                                conds = condsList[[2]], 
#                                target = targetList[[2]])

# hm1<-pheatmap(r2matlist3m$SeqPC2$r2mat_meta,cluster_cols = F,
#               cluster_rows=FALSE, fontsize = 8, 
#               silent=T,display_numbers = TRUE,
#               fontsize_number = 8,cellwidth = 20, cellheight = 20)
# hm2<-pheatmap(t(r2matlist3m$SeqPC2$r2mat_seq),
#               clustering_method="average",fontsize = 8,
#               silent = T,cluster_cols=FALSE, 
#               display_numbers = TRUE,fontsize_number = 8, cellwidth = 40)
# hms<-arrangeGrob(grobs = list(hm1[[4]],hm2[[4]]), 
#                  nrow = 1,widths = c(3,6))
# ggsave(plot = hms,filename = paste0(figdir,"hmsCorr1.pdf"),
#        width = 8,height = 5)


# # 3m cor1 - matrix --------------------
# pdf(paste0(figdir,"CorrPlot_Corr1N.pdf"),width = 9, height=8)
# pairs(cbind(PCAdatlist3m$SeqPC2$TopPC1,pairsdatList[[2]]),
#       col= condList[[2]],pch=19,upper.panel = panel.cor, 
#       cex.labels=1,gap = .25)
# dev.off()

# # 3m - cor1 - dea---------
# fitcontTisTime<-list()
# cornames<-paste0(tgrps,c("SeqPC4","SeqPC2"))
# datcorlist<-list(dat.cor1m$SeqPC4,
#                  dat.cor3m$SeqPC2)
# for(a in 1:2){
#   eset<-esetList[[a]]
#   an<-fData(eset)
#   f<-factor(eset$GT, levels=c("Hem","WT"))
#   design<-model.matrix(~0+f)
#   colnames(design)<-levels(f)
#   contrasts<-makeContrasts(Hem-WT, levels=design)
  
#   # uncorrected
#   datexpr<-normExprList[[a]]
#   fit <- lmFit(datexpr, design)
#   an<-an[rownames(datexpr),]
#   an$Gene.description<-gsub("[[:space:]]\\[.+","",an$Gene.description)
#   names(an)[names(an)=="Chromosome.scaffold.name"]<-"Chr"
#   fit$genes<-an[c("Gene.name","Gene.description")]
#   fit.cont <- contrasts.fit(fit, contrasts)
#   fit.cont <- eBayes(fit.cont)
#   rawname<-paste0(tgrps[a],"Pre")
#   fitcontTisTime[[rawname]]<-fit.cont
  
#   # corrected
#   datexpr<-datcorlist[[a]]
#   fit <- lmFit(datexpr, design)
#   an<-an[rownames(datexpr),]
#   an$Gene.description<-gsub("[[:space:]]\\[.+","",an$Gene.description)
#   names(an)[names(an)=="Chromosome.scaffold.name"]<-"Chr"
#   fit$genes<-an[c("Gene.name","Gene.description")]
#   fit.cont <- contrasts.fit(fit, contrasts)
#   fit.cont <- eBayes(fit.cont)
#   fitcontTisTime[[cornames[a]]]<-fit.cont
# }


# dtlistTisTime<-list()

# for(tgrp in names(fitcontTisTime)){
#   for(lfc in c(0,0.5,1,1.5)){
#     for(pval in c(0.05,0.1)){
#       dts<-as.data.frame(summary(decideTests(fitcontTisTime[[tgrp]],lfc = lfc,p.value = pval)))
#       dts$lfc<-lfc
#       dts$fdr<-pval
#       dts$grp<-tgrp
#       dtlistTisTime[[paste0(tgrp,"lfc",lfc,"p",pval)]]<-dts
#     }
#   }
  
# }


# dtdatTisTime<-do.call(rbind,dtlistTisTime)
# names(dtdatTisTime)[1:2]<-c("Dir","Cont")

# dtdatTisTimeCor<-dtdatTisTime[!dtdatTisTime$grp%in%c("str1mPre","str3mPre"),]
# dtdatTisTimeCor$lfc<-factor(dtdatTisTimeCor$lfc)
# ggplot(dtdatTisTimeCor,aes(lfc,Freq,fill=Dir))+
#   geom_col()+
#   facet_wrap(~grp,scales = "free")+theme(strip.text = element_text(size = 16),axis.text = element_text(size = 12))

# ggsave("dea_str.pdf",width = 8,height = 5)





# # dea------
# datexpr<-datExprfit
# GTtime<-factor(target$GTtime)
# design<-model.matrix(~0+GTtime)
# colnames(design)<-levels(GTtime)
# head(design)
# contrasts<-makeContrasts(Hem_1 - WT_1,Hem_3 - WT_3, levels=design)
# contrasts
# fit <- lmFit(datexpr, design)
# an<-fData(eset)
# an<-an[rownames(datexpr),]
# an$Gene.description<-gsub("[[:space:]]\\[.+","",an$Gene.description)
# names(an)[names(an)=="Chromosome.scaffold.name"]<-"Chr"
# fit$genes<-an[c("Gene.name","Gene.description")]
# fit.cont<-contrasts.fit(fit, contrasts)
# fit.cont<-eBayes(fit.cont)

# summary(decideTests(fit.cont))
# format(topTable(fit.cont,coef = "Hem_1 - WT_1" ,adjust="BH",p.value = .1, lfc=log2(1.5),number = 50),digits = 2)
# format(topTable(fit.cont,coef = "Hem_3 - WT_3" ,adjust="BH", p.value = .1, lfc=log2(1.5),number = 50),digits = 2)




