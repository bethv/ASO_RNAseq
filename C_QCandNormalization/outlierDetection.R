
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

dc2<-esets$dc[,!esets$dc$Mouse.ID%in%c("7648","7675")]
dc2<-countsGr5in30perfilt(dc2)
dc2<-exprs(dc2)
dc2<-log2(dc2+1)
x2<-as.data.frame(t(dc2))
rownames(x2)<-colnames(dc2)
dc2<-x2

dc.sampleTree = hclust(dist(datExprList$dc), method = "average")

dc2.sampleTree = hclust(dist(dc2), method = "average")

pdf(file = paste("dc","SampleClustering.pdf",sep = "_"),height = 5,width = 8)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(dc.sampleTree, main = paste("Sample clustering to detect outliers","dc", sep = " - "), sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

plot(dc2.sampleTree, main = paste("Sample clustering to detect outliers","dc2", sep = " - "), sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()

esets2<-list()
esets2$dc1m<-eset[,eset$tissue=="dc"&eset$Time==1]
esets2$dc3m<-eset[,eset$tissue=="dc"&eset$Time==3]
esets2$str1m<-eset[,eset$tissue=="str"&eset$Time==1]
esets2$str3m<-eset[,eset$tissue=="str"&eset$Time==3]


filtsets2<-lapply(esets2, countsGr5in30perfilt)

datExprList2<-lapply(filtsets2,function(x){
  x<-exprs(x)
  x<-log2(x+1)
  x2<-as.data.frame(t(x))
  rownames(x2)<-colnames(x)
  return(x2)
})

samptreeList<-lapply(datExprList2,function(x){
  hclust(dist(x), method = "average")
})


for (i in names(samptreeList)) {
  pdf(file = paste(i,"SampleClustering.pdf",sep = "_"),height = 5,width = 8)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(samptreeList[[i]], main = paste("Sample clustering to detect outliers",i, sep = " - "), sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
  dev.off()
}
