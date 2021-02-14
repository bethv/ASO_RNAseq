ff <- list.files( path = "./readcounts", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 )
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 2,drop=F] ) )
colnames(counts) <- gsub( "_ReadsPerGene\\.out\\.tab", "", gsub( "\\./readcounts/", "", ff ) )
row.names(counts) <-counts.files[[1]]$V1
saveRDS(counts, file="../data/star_counts.rds")


countinfo.files <- lapply( ff, read.table, nrows = 4 )
countinfo <- as.data.frame( sapply( countinfo.files, function(x) x[ , 2,drop=F] ) )
colnames(countinfo) <- gsub( "_ReadsPerGene\\.out\\.tab", "", gsub( "\\./readcounts/", "", ff ) )
row.names(countinfo) <-countinfo.files[[1]]$V1
saveRDS(countinfo, file="../data/star_countinfo.rds")
write.csv(countinfo,"../data/star_countinfo.csv" )
dim(counts)
counts<-counts[rowSums(counts)>0,]
row.names(counts)<-gsub("\\..+","",row.names(counts))
dim(counts)

pheno<-readRDS("../data/colDat.rds")
counts<-counts[,rownames(pheno)]
an<-readRDS("../data/mm_an.rds")
an<-an[rownames(counts),]

library(Biobase)
phenoData<-new("AnnotatedDataFrame", data=pheno) 
esetAll<- new("ExpressionSet", exprs=as.matrix(counts), phenoData=phenoData)
fData(esetAll)<-an

saveRDS(esetAll, file = "../../data/starCounts_esetAllSumGr0.rds")

