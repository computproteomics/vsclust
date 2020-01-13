########### functions for VSClust

validate <- shiny::validate

# extend to 702 cases:
LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)
cvalidate.xiebeni <- function(clres,m) {                                                
  xrows <- dim(clres$me)[1]                                                
  minimum <- -1                                                            
  error <- clres$within                                                    
  ncenters <- dim(clres$centers)[1]                                        
  for (i in 1:(ncenters - 1)) {                                            
    for (j in (i + 1):ncenters) {                                        
      diff <- clres$ce[i, ] - clres$ce[j, ]                            
      diffdist <- t(diff) %*% t(t(diff))                               
      if (minimum == -1)                                               
        minimum <- diffdist                                            
      if (diffdist < minimum)                                          
        minimum <- diffdist                                            
    }                                                                    
  }                                                                        
  xiebeni <- error/(xrows * minimum)                                       
  return(xiebeni)                                                          
}

# switch cluster numbers from largest to smallest
SwitchOrder <- function(Bestcl,NClust) {
  switching <- as.numeric(names(sort(table(Bestcl$cluster[rowMaxs(Bestcl$membership)>0.5]),decreasing=T)))
  if (length(switching)<NClust) 
    switching <- c(switching,which(!((1:NClust) %in% switching)))
  switching2 <- 1:NClust
  names(switching2) <- switching
  tBest <- Bestcl
  tBest$centers <- Bestcl$centers[switching,]
  rownames(tBest$centers) <- 1:NClust
  tBest$size <- Bestcl$size[switching]
  tBest$cluster <- switching2[as.character(Bestcl$cluster)]
  names(tBest$cluster) <- names(Bestcl$cluster)
  tBest$membership <- Bestcl$membership[,switching]
  colnames(tBest$membership) <- 1:NClust
  tBest
}


ClustComp <- function(tData,NSs=10,NClust=NClust,Sds=Sds, cores=1) {
  D<-ncol(tData)
  d<-sqrt(D/2)
  dims<-dim(tData)                                                                                                           
  mm <- 1 + (1418/dims[1] + 22.05)* dims[2]^(-2) + (12.33/dims[1] + 0.243)*dims[2]^(-0.0406*log(dims[1])-0.1134)
  
  ### d_i and d_t
  difunc <- function(c,D) { x <- 0:c; sum(choose(c,x)/(x*D+1)*(-1)^x) }
  
  di <- difunc(NClust,D)  /sqrt(pi) * gamma(D/2+1)^(1/D)
  dt <- (NClust)^(-1/D)
  
  p <- dnorm(di,0,Sds) * (1-dnorm(dt,0,Sds))^(NClust-1)
  
  m <- mm + p*mm*(D/3-1)
  m[m==Inf]<-0
  m[m==0]<-NA
  m[is.na(m)]<-mm*10
  #m<-rowMaxs(cbind(m,mm,na.rm=T))
  
  ## If m for highest Sd is mm then all = mm
  if (m[which.max(Sds)]== mm) 
    m[1:length(m)] <- mm
  
  #   plot(Sds,m,cex=0.5,pch=15,col=rainbow(maxClust)[NClust])
  #   abline(h=mm)
  
  colnames(tData)<-NULL
  PExpr <- new("ExpressionSet",expr=as.matrix(tData))
  PExpr.r <- filter.NA(PExpr, thres = 0.25)
  PExpr <- fill.NA(PExpr.r,mode = "mean")
  tmp <- filter.std(PExpr,min.std=0,visu=F)
  PExpr2 <- standardise(PExpr)
  
  cl <- makeCluster(cores)
  clusterExport(cl=cl,varlist=c("PExpr2","NClust","m"),envir=environment())
  clusterEvalQ(cl=cl, library(e1071FuzzVec))  
  clusterEvalQ(cl=cl, library(Biobase))  
  cls <- parLapply(cl,1:NSs, function(x) e1071FuzzVec::cmeans(exprs(PExpr2),NClust,m=m,verbose=F,iter.max=1000))
  print(cls[[1]])
  Bestcl <- cls[[which.min(lapply(cls,function(x) x$withinerror))]]
  cls <- parLapply(cl,1:NSs, function(x) e1071FuzzVec::cmeans(exprs(PExpr2),NClust,m=mm,verbose=F,iter.max=1000))
  Bestcl2 <- cls[[which.min(lapply(cls,function(x) x$withinerror))]]
  stopCluster(cl)
  
  # return validation indices
  list(indices=c(min(dist(Bestcl$centers)),cvalidate.xiebeni(Bestcl,mm),
                 min(dist(Bestcl2$centers)),cvalidate.xiebeni(Bestcl2,mm)),
       Bestcl=Bestcl,Bestcl2=Bestcl2,m=m,withinerror=Bestcl$withinerror,withinerror2=Bestcl2$withinerror) 
}


SignAnalPaired <- function(Data,NumCond,NumReps) {
  ##########################################################
  # significance analysis
  MAData<-Data[,2:(NumCond)]-Data[,1]
  for (i in 1:(NumReps-1))
    MAData<-cbind(MAData,Data[,(i*NumCond+1)+1:(NumCond-1)]-Data[,(i*NumCond+1)])
  rownames(MAData)<-rownames(Data)
  MAReps<-rep(1:(NumCond-1),NumReps)
  ##limma with ratios
  design<-plvalues<-NULL
  for (c in (1:(NumCond-1))) {
    design<-cbind(design,as.numeric(MAReps==c))
  }
  lm.fittedMA <- lmFit(MAData,design)
  lm.bayesMA<-eBayes(lm.fittedMA)
  topTable(lm.bayesMA)
  plvalues <- lm.bayesMA$p.value
  qvalues <- matrix(NA,nrow=nrow(plvalues),ncol=ncol(plvalues),dimnames=dimnames(plvalues))
  # qvalue correction
  for (i in 1:ncol(plvalues)) {
    tqs <- tryCatch(qvalue(na.omit(plvalues[,i]))$qvalues, 
                    error = function(e) NULL)
    if (length(tqs) >0) {
      qvalues[names(tqs),i] <- tqs
    }
    else {
      qvalues[names(tqs),i] <- NA
    }
  }
  
  return(list(plvalues=qvalues,Sds=sqrt(lm.bayesMA$s2.post)))
}

SignAnal <- function(Data,NumCond,NumReps) {
  ##########################################################
  # significance analysis
  Reps <- rep(1:NumCond,NumReps)
  design <- model.matrix(~0+factor(Reps-1))
  colnames(design)<-paste("i",c(1:NumCond),sep="")
  contrasts<-NULL
  First <- 1
  for (i in (1:NumCond)[-First]) contrasts<-append(contrasts,paste(colnames(design)[i],"-",colnames(design)[First],sep=""))
  contrast.matrix<-makeContrasts(contrasts=contrasts,levels=design)
  print(dim(Data))
  lm.fitted <- lmFit(Data,design)
  
  lm.contr <- contrasts.fit(lm.fitted,contrast.matrix)
  lm.bayes<-eBayes(lm.contr)
  topTable(lm.bayes)
  plvalues <- lm.bayes$p.value
  qvalues <- matrix(NA,nrow=nrow(plvalues),ncol=ncol(plvalues),dimnames=dimnames(plvalues))
  # qvalue correction
  for (i in 1:ncol(plvalues)) {
    tqs <- tryCatch(qvalue(na.omit(plvalues[,i]))$qvalues, 
                    error = function(e) NULL)
    print(tqs)
    if (length(tqs) >0) {
      qvalues[names(tqs),i] <- tqs
    }
    else {
      qvalues[names(tqs),i] <- NA
    }
    
  }
  return(list(plvalues=qvalues,Sds=sqrt(lm.bayes$s2.post)))
}


## overwrite function to set timeout limit higher
## OBSOLETE?
enrichDAVID2 <- function (gene, idType = "ENTREZ_GENE_ID", 
                          minGSSize = 5, maxGSSize = 500, annotation = "GOTERM_BP_ALL", pvalueCutoff = 0.05, 
                          pAdjustMethod = "BH", qvalueCutoff = 0.2, species = NA, david.user = "veits@bmb.sdu.dk") 
{
  Count <- List.Total <- Pop.Hits <- Pop.Total <- NULL
  pAdjustMethod <- match.arg(pAdjustMethod, c("bonferroni", 
                                              "BH"))
  david <- DAVIDWebService$new(email = david.user,url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  idType <- match.arg(idType, getIdTypes(david))
  setTimeOut(david,10000000)
  # browser()
  david.res <- addList(david, gene, idType = idType, listName = "clusterProfiler")
  if (david.res$inDavid == 0) {
    stop("All id can not be mapped. Please check 'idType' parameter...")
  }
  setAnnotationCategories(david, annotation)
  x <- getFunctionalAnnotationChart(david, threshold = 1, count = minGSSize)
  if (length(x@.Data) == 0) {
    warning("No significant enrichment found...")
    return(NULL)
  }
  term <- x$Term
  if (length(grep("~", term[1])) == 0) {
    sep <- ":"
  }
  else {
    sep <- "~"
  }
  term.list <- sapply(term, function(y) strsplit(y, split = sep))
  term.df <- do.call("rbind", term.list)
  # print(term.df)
  ID <- term.df[, 1]
  if (ncol(term.df)> 1) {
    Description <- term.df[, 2]
  } else {
    Description <- term.df[,1]
  }
  GeneRatio <- with(x, paste(Count, List.Total, sep = "/"))
  BgRatio <- with(x, paste(Pop.Hits, Pop.Total, sep = "/"))
  Over <- data.frame(ID = ID, Description = Description, GeneRatio = GeneRatio, 
                     BgRatio = BgRatio, pvalue = x$PValue)
  print(sum(duplicated(ID)))
  #row.names(Over) <- names(ID)
  if (pAdjustMethod == "bonferroni") {
    Over$p.adjust <- x$Bonferroni
  }
  else {
    Over$p.adjust <- x$Benjamini
  }
  qobj <- tryCatch(qvalue(p = Over$pvalue, lambda = 0.05, pi0.method = "bootstrap"), 
                   error = function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  }
  else {
    qvalues <- NA
  }
  Over$qvalue <- qvalues
  Over$geneID <- gsub(",\\s*", "/", x$Genes)
  Over$Count <- x$Count
  Over <- Over[Over$pvalue <= pvalueCutoff, ]
  Over <- Over[Over$p.adjust <= pvalueCutoff, ]
  if (!any(is.na(Over$qvalue))) {
    Over <- Over[Over$qvalue <= qvalueCutoff, ]
  }
  org <- getSpecieNames(david)
  org <- gsub("\\(.*\\)", "", org)
  # gc <- strsplit(Over$geneID, "/")
  if (!is.na(maxGSSize) || !is.null(maxGSSize)) {
    idx <- as.numeric(sub("/\\d+", "", Over$BgRatio)) <= 
      maxGSSize
    Over <- Over[idx, ]
  }
  
  # names(gc) <- Over$ID
  new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff, 
      pAdjustMethod = pAdjustMethod, organism = org, ontology = annotation, 
      gene = as.character(gene), keytype = idType)
}

#   calcBHI <- function(Accs,gos) {
#     ## enrichment does not yield all GO terms! This could lead to problems
#     BHI <- sumcomb <- vector("integer",length(Accs))
#     names(BHI) <- names(Accs)
#     for (i in names(gos@geneClusters)) {
#       genes <- Accs[[i]]
#       combs <- combn(genes,2)
#       sumcomb[i] <- ncol(combs)
#       clgroup <- gos@compareClusterResult[gos@compareClusterResult$Cluster==i,"geneID"]
#       for(j in 1:ncol(combs)) {
#         genepair <- combs[,j]
#         if(length(grep(genepair[2],grep(genepair[2],clgroup,value=T)))>0)
#           BHI[i] <- BHI[i] + 1
#       }
#     }
#     sum(BHI)/sum(sumcomb)
#   }


calcBHI <- function(Accs,gos) {
  ## enrichment does not yield all GO terms! This could lead to problems
  #  Rprof(tmp <- tempfile())
  BHI <- sumcomb <- vector("integer",length(Accs))
  names(BHI) <- names(Accs)
  goData <- as.data.frame(gos@compareClusterResult)
  for (i in names(gos@geneClusters)) {
    genes <- Accs[[i]]
    ispair <- matrix(F,length(genes),length(genes),dimnames = list(rows=genes,cols=genes))
    sumcomb[i] <- choose(length(genes),2)
    clgroup <- goData[goData$Cluster==i,"geneID"]
    for (j in 1:length(clgroup)) {
      tgenes <- genes[na.omit(match(unlist(strsplit(clgroup[j],"/"),use.names=F),as.character(genes)))]
      # print(tgenes)
      ltgenes <- length(tgenes)
      if (ltgenes > 1) {
        for (i1 in tgenes[1:(ltgenes-1)]) {
          # print(i1)
          ttgene <- tgenes[(which(i1==tgenes)+1):ltgenes]
          # for (i2 in tgenes[ttgene:ltgenes]) {
          # print(paste(tgenes[i1],tgenes[i2]))
          # ispair[i1,i2] <- ispair[i2,i1] <- T
          ispair[i1,ttgene] <- ispair[ttgene,i1] <- T
          # }
        }
        # combs <- combn(tgenes,2)
        
      }
      # print(tgenes)
    }
    BHI[i] <- BHI[i] + sum(ispair)/2
  }
  # Rprof()
  # print(summaryRprof(tmp))
  sum(BHI)/sum(sumcomb)
}


### Wrapper functions

## Wrapper for statistics
statWrapper <- function(dat, NumReps, NumCond, isPaired=F, isStat) {
  qvals <- statFileOut <- Sds <- NULL
  fdat <- dat
  tdat <- NULL
  print(isStat)
  if (isStat) {
    if(ncol(dat)!=NumReps*NumCond)
      print("Number of data columns must correspond to product of conditions and replicates!")
    if (isPaired) {
      print("Paired")
      ttt <- SignAnalPaired(dat,NumCond, NumReps)
    } else {
      ttt <- SignAnal(dat,NumCond, NumReps)
    }
    
    tdat<-rowMeans(dat[,seq(1,NumReps*NumCond,NumCond)],na.rm=T)
    
    print(Sds)
    
    Sds <- ttt$Sds
    qvals <- ttt$plvalues
    colnames(qvals) <- paste("qvalue ",LETTERS702[2:(NumCond)],"vsA",sep="")
    for (i in 2:NumCond) {
      tdat<-cbind(tdat,rowMeans(dat[,seq(i,NumReps*NumCond,NumCond)],na.rm=T))
    }
    colnames(tdat)<-paste("Mean of log ",LETTERS702[1:(NumCond)],sep="")
    dat <- cbind(tdat,Sds=Sds)
    #      dat <- dat[rowSums(is.na(ttt$plvalues))==0,]
    tdat <- dat[,1:(ncol(dat)-1)]
    Sds <- dat[,ncol(dat)]
  } else {
    Sds <- dat[,ncol(dat)]
    tdat <- dat[,1:(ncol(dat)-1)]
    NumReps <- 1
    NumCond <- ncol(dat)-1
  }
  dat[,ncol(dat)] <- dat[,ncol(dat)] / rowSds(as.matrix(tdat),na.rm=T)
  
  if (isStat) {
    statFileOut <- cbind(dat[,1:(ncol(dat)-1)],Sds,qvals)
  } else {
    statFileOut <- cbind(dat[,1:(ncol(dat)-1)],Sds)
  }
  
  # Remove columns with only 5% of the data
  plab <- rep(1:NumCond,NumReps)
  plab <- plab[colSums(!is.na(fdat))>0.05*nrow(fdat)]
  pcaDat <- fdat[,colSums(!is.na(fdat))>0.05*nrow(fdat)]
  # pcaDat <- t(pcaDat[complete.cases(pcaDat),])
  pcaDat <- (pcaDat[complete.cases(pcaDat),])
  validate(need(length(pcaDat)> 0, "Principal component analysis not shown as too many missing values"))      
  validate(need(nrow(pcaDat)> 10, "Principal component analysis not shown as too many missing values"))      
  pca<-prcomp(pcaDat,scale=T,retx=T)
  # scores <- pca$x
  # loadings <- pca$rotation
  scores <- pca$rotation
  loadings <- pca$x
  par(mfrow=c(1,2))
  tSds <- ifelse(is.na(dat[,ncol(dat)]),max(dat[,ncol(dat)],na.rm=T),sqrt(dat[,ncol(dat)]))
  print(max(as.integer(255/max(1/tSds)/tSds)))
  plot(loadings,cex=tSds,pch=16,col=paste("#000000",sprintf("%02X",as.integer(255/max(1/tSds)/tSds)),sep=""))
  title(main="Principal component analysis of data set (loadings)",sub="The point size corresponds to the estimated standard deviation")
  plot(scores,pch=19,col=rainbow(NumCond)[plab])
  title(main="Principal component analysis of data set (scores)",sub="Colors denote different conditions")
  legend("topright",paste("Condition",1:NumCond),pch=rep(19,NumCond),col=rainbow(NumCond)[1:NumCond])
  
  ## Preparing output
  Out <- list(dat=dat, qvals=qvals, statFileOut=statFileOut)
  Out
  
}

## Wrapper for estimation of cluster number
estimClustNum<- function(dat, maxClust=25, cores=1) {
  # print(head(rowSds(as.matrix(dat[,1:(ncol(dat)-1)]),na.rm=T)))
  Sds <- dat[,ncol(dat)] / rowSds(as.matrix(dat[,1:(ncol(dat)-1)]),na.rm=T)
  ClustInd<-matrix(NA,nrow=maxClust,ncol=6)
  multiOut <- lapply(3:maxClust,function(x) {    
    if (!is.null(getDefaultReactiveDomain())) {
      incProgress(1, detail = paste("Running cluster number",x))
    } else {
      print(paste("Running cluster number",x))
    }
    clustout <- ClustComp(dat[,1:(ncol(dat)-1)],NClust=x,Sds=dat[,ncol(dat)],NSs=16, cores)
    c(clustout$indices,sum(rowMaxs(clustout$Bestcl$membership)>0.5),
      sum(rowMaxs(clustout$Bestcl2$membership)>0.5))
  })
  for (NClust in 3:maxClust) 
    ClustInd[NClust,] <- multiOut[[NClust-2]]
  # for (NClust in 3:maxClust) {
  #   print(paste("Running cluster number", NClust))
  #   if (!is.null(getDefaultReactiveDomain())) {
  #     incProgress(1, detail = paste("Running cluster number",NClust))
  #   } else {
  #     print(paste("Running cluster number",NClust))
  #   }
  #   clustout <- ClustComp(dat[,1:(ncol(dat)-1)],NClust=NClust,Sds=dat[,ncol(dat)],NSs=16, cores)
  #   ClustInd[NClust,]<-c(clustout$indices,sum(rowMaxs(clustout$Bestcl$membership)>0.5),
  #                        sum(rowMaxs(clustout$Bestcl2$membership)>0.5))
  # }
  
  dmindist <- c(which.max(ClustInd[3:(maxClust-2),1]-ClustInd[4:(maxClust-1),1])+2,
                which.max(ClustInd[3:(maxClust-2),3]-ClustInd[4:(maxClust-1),3])+2)
  dxiebeni <- c(which.min(ClustInd[3:(maxClust-1),2])+2,
                which.min(ClustInd[3:(maxClust-1),4])+2)  
  print(dmindist)
  plot.new() ## clean up device
  par(mfrow=c(1,3))
  plot(3:maxClust,ClustInd[3:(maxClust),1],col=2 , type="b", 
       main="Min. centroid distance\n(Highest jump is best)",
       xlab="Number of clusters", ylab="Index",ylim=c(min(ClustInd[,c(1,3)],na.rm=T),max(ClustInd[,c(1,3)],na.rm=T)))
  lines(3:maxClust,ClustInd[3:(maxClust),3],col=3,type="b")
  points(dmindist[1],ClustInd[dmindist[1],1],pch=15,col=1,cex=2)
  legend("topright",legend = c("VSClust","Standard"), lty=c(1,1),col=2:3)
  grid(NULL,NA,lwd=1,col=1)
  plot(3:maxClust,ClustInd[3:(maxClust),2], col=2, type="b", main="Xie-Beni index\n(Lowest is best)",
       xlab="Number of clusters", ylab="Index",ylim=c(min(ClustInd[,c(2,4)],na.rm=T),max(ClustInd[,c(2,4)],na.rm=T)))
  lines(3:maxClust,ClustInd[3:(maxClust),4],type="b",col=3)
  points(dxiebeni[1],ClustInd[dxiebeni[1],2],pch=15,col=1,cex=2)
  legend("topright",legend = c("VSClust","Standard"), lty=c(1,1),col=2:3)
  grid(NULL,NA,lwd=1,col=1)
  plot(3:maxClust,ClustInd[3:(maxClust),5], col=2, type="b", main="Total number of assigned features",
       xlab="Number of clusters", ylab="Assigned features",ylim=c(min(ClustInd[,5:6],na.rm=T),max(ClustInd[,5:6],na.rm=T)))
  lines(3:maxClust,ClustInd[3:(maxClust),6],type="b",col=3)
  legend("topright",legend = c("VSClust","Standard"), lty=c(1,1),col=2:3)
  # finally plot
  p <- recordPlot()
  
  # Output
  Out <- list(ClustdInd=ClustInd, p=p, numclust=dmindist[1])
  Out
}

## Wrapper for clustering
runClustWrapper <- function(dat, NClust, proteins=NULL, VSClust=T, cores) {
  # dat <- dat[rowSums(is.na(dat))==0,]
  PExpr <- new("ExpressionSet",expr=as.matrix(dat[,1:(ncol(dat)-1)]))
  PExpr.r <- filter.NA(PExpr, thres = 0.25)
  PExpr <- fill.NA(PExpr.r,mode = "mean")
  tmp <- filter.std(PExpr,min.std=0,visu=F)
  PExpr <- standardise(PExpr)
  
  
  clustout <- ClustComp(exprs(PExpr),NClust=NClust,Sds=dat[,ncol(dat)],NSs=16, cores)
  if (VSClust) {
    Bestcl <- clustout$Bestcl
  } else {
    Bestcl <- clustout$Bestcl2
  }
  Bestcl <- SwitchOrder(Bestcl,NClust)
  
  # sorting for membership values (globally)
  Bestcl$cluster <- Bestcl$cluster[order(rowMaxs(Bestcl$membership,na.rm=T))]
  Bestcl$membership <- Bestcl$membership[order(rowMaxs(Bestcl$membership,na.rm=T)),]
  PExpr <- PExpr[names(Bestcl$cluster),]
  
  if (!is.null(getDefaultReactiveDomain()))
    incProgress(0.7, detail = paste("Plotting",NClust))
  
  plot.new() ## clean up device
  par(lwd=0.25)
  oldmar <- par("mar")
  par(mar=c(2,2,3,3),mgp=c(2,1,0))
  par(mar=par("mar")/max(1,NClust/20))
  mfuzz.plot2(PExpr,cl=Bestcl,mfrow=c(round(sqrt(NClust)),ceiling(sqrt(NClust))),min.mem=0.5,x11=F,colo="fancy")
  # Mfuzz::mfuzzColorBar(col="fancy")
  p <- recordPlot()
  par(lwd=1,mar=oldmar)
  
  colnames(Bestcl$membership) <- paste("membership of cluster",colnames(Bestcl$membership))
  outFileClust <- exprs(PExpr)
  if (!is.null(proteins)) {
    outFileClust <- cbind(outFileClust,names=as.character(proteins[rownames(outFileClust)]))
  }
  
  rownames(Bestcl$centers) <- paste("Cluster",rownames(Bestcl$centers))
  ClustInd <- as.data.frame(table(Bestcl$cluster[rowMaxs(Bestcl$membership)>0.5]))
  if (ncol(ClustInd) == 2)
    colnames(ClustInd) <- c("Cluster","Members")
  else 
    ClustInd <- cbind(1:max(Bestcl$cluster),rep(0,max(Bestcl$cluster)))
  
  ## Output
  Out <- list(dat=PExpr, Bestcl=Bestcl, p=p, outFileClust=outFileClust, ClustInd=ClustInd)
  return(Out)
}

# Wrapper for functional enrichment (TODO)
runFuncEnrich <- function(cl, dat, protnames, idtypes, infosource) {
  Accs <- list()
  for (c in 1:max(cl$cluster)) {
    Accs[[c]] <- names(which(cl$cluster==c & rowMaxs(cl$membership)>0.5))
    Accs[[c]] <- Accs[[c]][Accs[[c]]!=""]
    if (length(Accs[[c]])>0) {
      
      if (length(Accs[[c]])>1) {
        tdat <- (dat[Accs[[c]],])
      } else {
        tdat <- t(dat[Accs[[c]],])
      }
      if (!is.null(protnames)) {
        Accs[[c]] <- as.character(protnames[Accs[[c]]])
      }
      
      Accs[[c]] <- sub("-[0-9]","",Accs[[c]])
    } else {
      tdat <- t(rep(NA,ncol(dat)+1))
    }
  }
  # TODO: add extraction of multiple accession numbers
  names(Accs) <- paste("Cluster",1:length(Accs))
  Accs <- lapply(Accs,function(x) unique(ifelse(is.na(x),"B3",x)))
  Accs <- Accs[lapply(Accs,length)>0]
  x <- NULL
  try(x <- compareCluster(Accs, fun="enrichDAVID", annotation=infosource,
                          idType=idtypes,
                          david.user = "veits@bmb.sdu.dk"))
  validate(need(!is.null(x),"No result. Wrong ID type?"))
  incProgress(0.7, detail = "received")
  print("got it")
  x@compareClusterResult <- cbind(x@compareClusterResult,log10padval=log10(x@compareClusterResult$p.adjust))
  print(x@compareClusterResult)
  y <- new("compareClusterResult",compareClusterResult=x@compareClusterResult)
  if (length(unique(y@compareClusterResult$ID)) > 20) {
    print("Reducing number of DAVID results")
    y@compareClusterResult <- y@compareClusterResult[
      order(y@compareClusterResult$p.adjust)[1:20],]
    y@compareClusterResult$Cluster <- as.character(y@compareClusterResult$Cluster)
    print(x)
  }
  BHI <- calcBHI(Accs,x)
  return(list(fullFuncs=x, redFuncs=y, BHI=BHI))
  
}