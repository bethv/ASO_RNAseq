


dlines<-function(counts,log=TRUE,main=""){
  if(log==TRUE){ 
    dat <- log(counts[,1],10)
    xlab="Raw read counts per gene (log10)"
  }else{
    dat<-counts[,1]
    xlab="Raw read counts per gene"
  }
  d <- density(dat)
  plot(d,xlim=c(1,8),main=main,ylim=c(0,.45),xlab=xlab, ylab="Density",cex.main=2)
  
  for (s in 2:dim(counts)[2]){
    if(log==TRUE){
      dat <- log(counts[,s],10) 
    }else{
      dat <- counts[,s]
    }
    d <- density(dat)
    lines(d)
  } 
}


dboxes<-function(counts, log=TRUE,main=""){
  if(log==TRUE){
    dat <- log(counts,10)
    ylab="Raw read counts per gene (log10)"
  }else{
    dat<-counts
    ylab="Raw read counts per gene"
  }
  suppressWarnings(boxplot(dat, main=main, xlab="", ylab=ylab,axes=FALSE))
  axis(2)
  axis(1,at=c(1:dim(counts)[2]),labels=colnames(counts),las=2,cex.axis=0.8)
  
}


getplotdata<-function(data,conds2,kfit,topn=NULL){
  cv <- as.matrix(conds2)
  dim(cv) <- c(1,prod(dim(cv)))
  names <- rownames(table(cv))
  col_default <- c("red","blue","green","orange","bisque4", "black",
                   "brown","cyan", "darkgreen", "darkgrey", "darkmagenta",
                   "darkolivegreen", "darkorange", "darkred", "darkslateblue",
                   "darkturquoise", "floralwhite", "greenyellow","grey", 
                   "lightcyan", "lightcyan1", "lightgreen", "lightsteelblue1",
                   "lightyellow","magenta", "mediumpurple3","midnightblue",
                   "paleturquoise", "pink", "plum1", "plum2", "royalblue",
                   "saddlebrown", "salmon", "sienna3", "skyblue", "skyblue3",
                   "steelblue", "tan", "thistle1", "thistle2", "turquoise",
                   "violet", "white", "yellowgreen", "grey60 ", "orangered4",
                   "brown4", "darkorange2", "ivory")
  col <- rainbow(length(names))
  names(col)<- names
  clab <- matrix(col[as.matrix(conds2)],nrow=nrow(conds2),ncol=ncol(conds2))
  colnames(clab) <- colnames(conds2) 
  rownames(clab)<-rownames(conds2)
  n=dim(conds2)[2]
  for(k in 1:n){
    a=1;for (i in unique(clab[,k])){
      clab[which(clab[,k]==i),k]=col_default[a]
      a=a+1
    }
  }
  if(!is.null(topn)){
    rowvar <- apply(data,1,var)
    ordV <- order(rowvar,decreasing=TRUE)
    topset <- head(ordV, n=topn)
    data <- data[topset,]
    main=paste0("Top ",topn, " by variance")
  }else{
    main="All genes"
  }
  ldat <- dist(t(data))
  invisible(fit <- isoMDS(ldat, k=kfit,))
  
  return(list(fit=fit,clab=clab, conds2=conds2, main=main))
}


makeplot<-function(plotdat,c1,c2,main=NULL,legpos="topleft",
                   legH=FALSE,legtxtw=0.045,legxint=0.25,
                   axes=T,frame=T,sizeText=get("sizeText",envir = parent.frame())){
  fit=plotdat$fit
  conds2=plotdat$conds2
  clab=plotdat$clab
  if(!is.null(main)){
    plot.title<-main
  }else{
    plot.title<-plotdat$main
  }
  x <- fit$points[,c1]; y <- fit$points[,c2]
  for(i in 2:dim(conds2)[2]){
    plot(x, y, xlab=paste("Coordinate",c1),
         ylab=paste("Coordinate",c2),
         main=plot.title, type="p", pch=16, cex=1,
         col=clab[,i], xlim=c(min(x)*1.4, max(x)*1.4),
         ylim=c(min(y)*1.4, max(y)*1.4),axes=axes, frame.plot=frame)
    text(x,y+(max(y)-min(y))/50, labels = conds2$sample, cex=sizeText,col=clab[,i])
    
    if(legH==TRUE){
      legend(legpos, legend=unique(as.factor(conds2[,i])),fill=unique(clab[,i]),cex=sizeText, 
             horiz=TRUE,x.intersp=legxint, text.width=legtxtw,bty = "n")
    }else{
      legend(legpos, legend=unique(as.factor(conds2[,i])),fill=unique(clab[,i]),cex=sizeText,bty = "n")
    }
    
  }
}





# makeplot<-function(c1,c2,fit,conds2,clab,main,legpos="topleft",legH=FALSE,legtxtw=0.045,legxint=0.25){
#   x <- fit$points[,c1]; y <- fit$points[,c2]
#   for(i in 2:dim(conds2)[2]){
#     plot(x, y, xlab=paste("Coordinate",c1),
#          ylab=paste("Coordinate",c2),
#          main=main, type="p", pch=16, cex=1,
#          col=clab[,i], xlim=c(min(x)*1.4, max(x)*1.4),
#          ylim=c(min(y)*1.4, max(y)*1.4))
#     text(x,y+(max(y)-min(y))/50, labels = conds2$sample, cex=sizeText,col=clab[,i])
#     
#     if(legH==TRUE){
#       legend(legpos, legend=unique(as.factor(conds2[,i])),fill=unique(clab[,i]),cex=sizeText, 
#              horiz=TRUE,x.intersp=legxint, text.width=legtxtw,bty = "n")
#     }else{
#       legend(legpos, legend=unique(as.factor(conds2[,i])),fill=unique(clab[,i]),cex=sizeText,bty = "n")
#     }
#     
#   }
# }

# function to get expression PCs
# PCvar<-function(data){
#   thisdat.HTSC <- t(scale(t(data),scale=F))
#   PC.HTSC <- prcomp(thisdat.HTSC,center=F);
#   TopPC1 <- PC.HTSC$rotation[,1:5];
#   varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
#   topvar <- varexp[1:10]
#   colnames(TopPC1) <- paste("Exp.PC\n",
#                             colnames(TopPC1)," (", signif(100*topvar[1:5],2),"%)",sep="")
#   return(list(TopPC1=TopPC1,
#               varexp=varexp,
#               topvar=topvar))
# }

#function to plot heatmaps
get.r2mats<-function(topPC,conds,target){
  n=dim(conds)[2]
  r2mat_meta = matrix(NA,nrow=5,ncol=n-1)
  r2mat_seq = matrix(NA,nrow=5,ncol=dim(target)[2]-n)
  rownames(r2mat_meta) <- rownames(r2mat_seq) <- colnames(topPC)
  datMeta_model = target[,2:n,drop=F]
  datSeq_model=target[,-c(1:n)]
  
  colnames(r2mat_meta) <- colnames(datMeta_model)
  colnames(r2mat_seq) <- colnames(datSeq_model)
  for(i in c(1:5)){
    for(j in c(1:dim(datMeta_model)[2])){
      to_remove = which(is.na(datMeta_model[,j])==TRUE)
      if(length(to_remove)>0){
        tmp_topPC=topPC[-to_remove,i]
        tmp_meta=datMeta_model[-to_remove,j]}
      else{
        tmp_topPC=topPC[,i]
        tmp_meta=datMeta_model[,j]  
      }
      mod_mat=model.matrix(~tmp_meta)[,-1]
      mod=summary(lm(tmp_topPC~mod_mat))
      r2mat_meta[i,j]=mod$adj.r.squared
    }
  }   
  
  for(i in c(1:5)){
    for(j in c(1:dim(datSeq_model)[2])){
      to_remove = which(is.na(datSeq_model[,j])==TRUE)
      if(length(to_remove)>0){
        tmp_topPC=topPC[-to_remove,i]
        tmp_seq=datSeq_model[-to_remove,j]}
      else{
        tmp_topPC=topPC[,i]
        tmp_seq=datSeq_model[,j]  
      }
      mod=summary(lm(tmp_topPC~as.numeric(tmp_seq)))
      r2mat_seq[i,j]=mod$adj.r.squared
    }
  }
  
  return(list(r2mat_meta=r2mat_meta,
              r2mat_seq=r2mat_seq,
              datMeta_model=datMeta_model,
              datSeq_model=datSeq_model))
  
}

# PCvar<-function(data){
#   thisdat.HTSC <- t(scale(t(data),scale=F))
#   PC.HTSC <- prcomp(thisdat.HTSC,center=F);
#   TopPC1 <- PC.HTSC$rotation[,1:5];
#   varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
#   topvar <- varexp[1:10]
#   colnames(TopPC1) <- paste("Exp.PC\n",
#                             colnames(TopPC1)," (", signif(100*topvar[1:5],2),"%)",sep="")
#   return(list(TopPC1=TopPC1,
#               varexp=varexp,
#               topvar=topvar))
# }

PCvar<-function(data){
  thisdat.HTSC <- t(scale(t(data),scale=F))
  PC.HTSC <- prcomp(thisdat.HTSC,center=F);
  TopPC1 <- PC.HTSC$rotation[,1:5];
  varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
  topvar <- varexp[1:10]
  colnames(TopPC1) <- paste("Exp_",
                            colnames(TopPC1),"\n(", signif(100*topvar[1:5],2),"%)",sep="")
  return(list(TopPC1=TopPC1,
              varexp=varexp,
              topvar=topvar))
}

#function to plot correlations
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (class(x) == "numeric" & class(y) == "numeric") {
    r <- abs(cor(x, y,use="pairwise.complete.obs",method="spearman"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$adj.r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


makemds.set<-function(data=get("normDat", envir = parent.frame()),
                      conds2=get("conds2", envir = parent.frame()),
                      kfit=6
                      ){
  mdlist<-list(All = getplotdata(data = data, conds2 = conds2, kfit = kfit,topn = NULL),
               top500 = getplotdata(data = data, conds2 = conds2, kfit = kfit,topn = 500),
               top1000 = getplotdata(data = data, conds2 = conds2, kfit = kfit,topn = 1000)
  )
  
  return(mdlist)
  
}

mdsplot1v2<-function(mds.set,
                     outdir=get("figdir", envir = parent.frame()),
                     fname){
  pdf(paste0(outdir,fname),width = 10, height = 5)
  sizeText=1; par(mfrow=c(1,3))
  for(i in 1:length(mds.set)){
    plotdat<-mds.set[[i]]
    makeplot(plotdat = plotdat,1,2,)
  }
  dev.off()
}

mdsplotAll<-function(mds.set,
                     outdir=get("figdir", envir = parent.frame()),
                     fname){
  pdf(paste0(outdir,fname),width = 10, height = 5)
  sizeText=.75; par(mar = c(0,0,0,0))
  layout(mat = matrix(seq(1,16),ncol = 4, byrow = T),
         heights = c(.5,4,4,4),widths = c(2,4,4,4))
  plot.new()
  plot.new();text(0.5,0.5,"All",cex=2)
  plot.new();text(0.5,0.5,"Top 500 (var)",cex=2)
  plot.new();text(0.5,0.5,"Top 1000 (var)",cex=2)
  plot.new();text(0.5,0.5,"1 v 2",cex=2)
  for(i in 1:length(mds.set)){
    plotdat<-mds.set[[i]]
    makeplot(plotdat = plotdat,1,2,main="",axes = F)
  }
  plot.new();text(0.5,0.5,"3 v 4",cex=2)
  for(i in 1:length(mds.set)){
    plotdat<-mds.set[[i]]
    makeplot(plotdat = plotdat,3,4,main="",axes = F)
  }
  plot.new();text(0.5,0.5,"5 v 6",cex=2)
  for(i in 1:length(mds.set)){
    plotdat<-mds.set[[i]]
    makeplot(plotdat = plotdat,5,6,main="",axes = F)
  }
  dev.off()
}


mdsCompPlot<-function(mds.list,
                      outdir=get("figdir", envir = parent.frame()),
                      fname, fw=10, fh=5){
  pdf(paste0(outdir,fname),width = fw, height = fh)
  sizeText=.75; par(mar = c(0,0,0,0))
  layout(mat = matrix(seq(1,(length(mds.list)+1)*4),ncol = 4,byrow = T),
         heights = c(.5,rep(4,length(mds.list))),widths = c(2,4,4,4))
  plot.new();text(0.5,0.5,"Coord 1 v 2",cex=2)
  plot.new();text(0.5,0.5,"All",cex=2)
  plot.new();text(0.5,0.5,"Top 500 (var)",cex=2)
  plot.new();text(0.5,0.5,"Top 1000 (var)",cex=2)
  
  for(i in 1:length(mds.list)){
    plot.new();text(0.5,0.5, names(mds.list)[i],cex=2)
    for(j in 1:length(mds.list[[i]])){
      plotdat<-mds.list[[i]][[j]]
      makeplot(plotdat = plotdat,1,2,main="",axes = F)
    }
  }
  dev.off()
}







plotvar<-function(pclist,expPCs,
                  outdir=get("figdir", envir = parent.frame()),
                  fname){
  
  pdf(paste0(outdir,fname),height = 4,width = 7.5)
  par(mfrow=c(1,2))
  plot(pclist[[expPCs]][["topvar"]], xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", 
       ylab = "variance explained", main="Variance explained by Expression PCs",cex.main=1)
  plot(pclist$seq$topvar, xlim = c(0, 10), type = "b", pch = 16, xlab = "principal components", 
       ylab = "variance explained", main="Variance explained by sequencing PCs",cex.main=1)
  dev.off()
}

plothms<-function(pclist,expPCs,
                  conds=get("conds", envir = parent.frame()),
                  target=get("target", envir = parent.frame()), 
                  outdir=get("figdir", envir = parent.frame()),
                  fname){
  r2mats<-get.r2mats(pclist[[expPCs]][["TopPC1"]],conds = conds,target = target)
  hm1<-pheatmap(r2mats$r2mat_meta,cluster_cols = F,cluster_rows=FALSE, fontsize = 8, 
                silent=T,display_numbers = TRUE,fontsize_number = 8,cellwidth = 20, cellheight = 20)
  hm2<-pheatmap(t(r2mats$r2mat_seq),clustering_method="average",fontsize = 8,silent = T,
                cluster_cols=FALSE, display_numbers = TRUE,fontsize_number = 8, cellwidth = 40)
  hms<-arrangeGrob(grobs = list(hm1[[4]],hm2[[4]]), 
                   nrow = 1,widths = c(3,6))
  ggsave(plot = hms,filename = paste0(outdir,fname),
         width = 8,height = 5)
  
}

cormatrix<-function(pclist,expPCs,
                    pairsdat=get("pairsdat", envir = parent.frame()),
                    cond=get("cond", envir = parent.frame()), 
                    outdir=get("figdir", envir = parent.frame()), fname, N=T){
  if(N==TRUE){
    fw = 9; fh=8; cex.lab=1
    fname<-gsub("\\.pdf","N\\.pdf",fname)
  }else{
    fw = 18; fh=8; cex.lab=1.5
  }
  pdf(paste0(outdir,fname),width = fw, height=fh)
  pairs(cbind(pclist[[expPCs]][["TopPC1"]],pairsdat),
        col= cond,pch=19,upper.panel = panel.cor, cex.labels=cex.lab,gap = .25)
  dev.off()
}


correct<-function(modVars=get("modVars", envir = parent.frame()), 
                  Y=get("normExpr",envir = parent.frame()),
                  regColIdx=ifelse("Time"%in%modVars,4,3),interaction=FALSE){
  model<-paste("~",paste(modVars, collapse="+"))
  print(model)
  X =eval(parse(text = paste0("model.matrix(",model,")")))
  endIndex<-ifelse(interaction==FALSE,dim(X)[2],dim(X)[2]-1)
  print(head(X))
  print(paste0("cols to correct for: ",regColIdx,":",endIndex))
  beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  b = as.data.frame(t(beta))
  to_regress = (as.matrix(X[,regColIdx:endIndex,drop=FALSE]) %*% 
                  (as.matrix(beta[regColIdx:endIndex,,drop=FALSE])))
  datExprfit = normExpr - t(to_regress)
  return(datExprfit)
}



