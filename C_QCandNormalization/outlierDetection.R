
library(Biobase)
library(WGCNA)
options(stringsAsFactors = FALSE)

eset<-readRDS("data/rawCountEset.rds")

esets<-list()
esets$dc<-eset[,eset$tissue=="dc"]
esets$str<-eset[,eset$tissue=="str"]

countsGr5in30perfilt<-function(eset){
  n30percent<-.3*dim(eset)[2]
  countsGr5in30per<-apply(exprs(eset),1,function(y){
    y2<-ifelse(y>5,1,0)
    y3<-ifelse(sum(y2)>n30percent,TRUE,FALSE)
    return(y3)
  })
  esetfilt<-eset[countsGr5in30per,]
  return(esetfilt)
}

filtsets<-lapply(esets, countsGr5in30perfilt)


datExprList<-lapply(filtsets,function(x){
  x<-exprs(x)
  x<-log2(x+1)
  x2<-as.data.frame(t(x))
  rownames(x2)<-colnames(x)
  return(x2)
})

for (i in names(datExprList)) {
  datExpr0<-datExprList[[i]]
  sampleTree = hclust(dist(datExpr0), method = "average")
  pdf(file = paste(i,"SampleClustering.pdf",sep = "_"),height = 5,width = 8)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = paste("Sample clustering to detect outliers",i, sep = " - "), sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
  dev.off()
}
