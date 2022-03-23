#' vsclust: variance-senstive clustering
#'
#' @section VSClust functions:
#' The mypackage functions ...
#'
#' @docType package
#' @name vsclust
#' @useDynLib vsclust
NULL
#> NULL

########### functions for VSClust


validate <- shiny::validate

# extend to 702 cases:
LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))


# Some mathematical functions
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)



#' Calculate the Xie Beni index for validation of the cluster number
#'
#' @param clres Output from clustering. Either fclust object or similar list
#' @param m Fuzzifier value
#' @return Xie Beni index
#' @examples
#' TODO
#' @export
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

#' Run the vsclust clustering algorithm
#'
#' @param x a numeric data matrix 
#' @param centers Either numeric for number of clusters or numeric matrix with center coordinates
#' @param iter.max Numeric for maximum number of iterations
#' @param verbose Verbose information
#' @param m Fuzzifier value: numeric or vector of length equal to number of rows of x 
#' @param rate.par (experimental) numeric value for punishing missing values
#' @param weights numeric or vector of length equal to number of rows of x 
#' @param control list with arguments to vsclust algorithms (now only cutoff for relative tolerance: reltol)
#' @return Xie Beni index
#' @examples
#' TODO
#' @export
#' @useDynLib
vsclust_algorithm <-
  function(x, centers, iter.max = 100, verbose = FALSE,
           dist = "euclidean", m = 2,
           rate.par = NULL, weights = 1, control = list())
  {
    x <- as.matrix(x)
    xrows <- nrow(x)
    xcols <- ncol(x)
    
    if(missing(centers))
      stop("Argument 'centers' must be a number or a matrix.")
    
    dist <- pmatch(dist, c("euclidean", "manhattan"))
    if(is.na(dist)) 
      stop("invalid distance")
    if(dist == -1) 
      stop("ambiguous distance")
    
    if(length(centers) == 1) {
      ncenters <- centers
      centers <- x[sample(1 : xrows, ncenters), , drop = FALSE]
      centers[is.na(centers)] <- 0
      if(any(duplicated(centers))) {
        cn <- unique(x)
        mm <- nrow(cn)
        if(mm < ncenters) 
          stop("More cluster centers than distinct data points.")
        centers <- cn[sample(1 : mm, ncenters), , drop = FALSE]
      }
    } else {
      centers <- as.matrix(centers)
      if(any(duplicated(centers))) 
        stop("Initial centers are not distinct.")
      cn <- NULL
      ncenters <- nrow(centers)
      if (xrows < ncenters)
        stop("More cluster centers than data points.")
    }
    
    if(xcols != ncol(centers))
      stop("Must have same number of columns in 'x' and 'centers'.")
    
    if(iter.max < 1) 
      stop("Argument 'iter.max' must be positive.")
    
    if(missing(rate.par)) {
      rate.par <- 0
    }
    
    reltol <- control$reltol
    if(is.null(reltol))
      reltol <- sqrt(.Machine$double.eps)
    if(reltol <= 0)
      stop("Control parameter 'reltol' must be positive.")
    
    if(any(weights < 0))
      stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
      stop("Argument 'weights' has no positive elements.")
    weights <- rep(weights, length = xrows)
    weights <- weights / sum(weights)
    
    # if length of fuzzifiers is lower than number of features, repeat the pattern until end. 
    # Also counts for single fuzzifier
    m <- rep(m, length = xrows)
    
    initcenters <- centers
    pos <- as.factor(1 : ncenters)
    rownames(centers) <- pos
    
    u <- matrix(0.0,  nrow = xrows, ncol = ncenters)
    iter <- c(0L)
    val <- vsclust:::c_plusplus_means(x, centers, weights, m, dist-1, iter.max, reltol, verbose, u , 1, iter, NA, rate.par)
    # put modified values in retval
    retval <- list(x=x, xrows = xrows, xcols = xcols, centers = centers,ncenters=ncenters, m = m, dist = dist -1,
                   iter.max = iter.max, reltol = reltol, verbose = verbose, rate.par = rate.par, u = u, ermin = val, iter = iter)
    
    centers <- matrix(retval$centers, ncol = xcols,
                      dimnames = list(1 : ncenters,
                                      colnames(initcenters)))
    u <- matrix(retval$u, ncol = ncenters,
                dimnames = list(rownames(x), 1 : ncenters))
    # u <- u[order(perm), ]
    iter <- retval$iter - 1
    withinerror <- retval$ermin
    
    cluster <- apply(u, 1, which.max)
    clustersize <- as.integer(table(cluster))
    
    retval <- list(centers = centers, size = clustersize,
                   cluster = cluster, membership = u, iter = iter,
                   withinerror = withinerror, call = match.call())
    class(retval) <- c("fclust")
    return(retval)
  }

# Function to return content of fclust object
print.fclust <-
  function(x, ...)
  {
    cat("Error:\n")
    print(x$withinerror, ...)
    cat("Fuzzy c-means clustering with", length(x$size), "clusters\n")
    cat("\nCluster centers:\n")
    print(x$centers, ...)
    cat("\nMemberships:\n")
    print(x$membership, ...)
    cat("\nClosest hard clustering:\n")
    print(x$cluster, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
  }

#' switch cluster numbers from largest to smallest
#' @importFrom matrixStats rowMaxs
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

#' Running VSClust on given data set (with data pre-processing)
#' @import parallel
#' @export
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
  
  ## If m for highest Sd is mm then all = mm
  if (m[which.max(Sds)]== mm) 
    m[1:length(m)] <- mm
  
  cl <- makeCluster(cores)
  clusterExport(cl=cl,varlist=c("tData","NClust","m","vsclust_algorithm"),envir=environment())
  clusterEvalQ(cl=cl, library(vsclust))  
  #clusterEvalQ(cl=cl, library(Biobase))  
  cls <- parLapply(cl,1:NSs, function(x) vsclust_algorithm(tData,NClust,m=m,verbose=F,iter.max=1000))
  #print(cls[[1]])
  Bestcl <- cls[[which.min(lapply(cls,function(x) x$withinerror))]]
  cls <- parLapply(cl,1:NSs, function(x) vsclust_algorithm(tData,NClust,m=mm,verbose=F,iter.max=1000))
  Bestcl2 <- cls[[which.min(lapply(cls,function(x) x$withinerror))]]
  stopCluster(cl)
  
  # return validation indices
  list(indices=c(min(dist(Bestcl$centers)),cvalidate.xiebeni(Bestcl,mm),
                 min(dist(Bestcl2$centers)),cvalidate.xiebeni(Bestcl2,mm)),
       Bestcl=Bestcl,Bestcl2=Bestcl2,m=m,withinerror=Bestcl$withinerror,withinerror2=Bestcl2$withinerror) 
}

#' Statistical analysis (paired)
#' @import limma
#' @importFrom qvalue qvalue
#' @export
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

#' Statistical analysis (paired)
#' @import limma
#' @importFrom qvalue qvalue
#' @export
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
  #topTable(lm.bayes)
  plvalues <- lm.bayes$p.value
  qvalues <- matrix(NA,nrow=nrow(plvalues),ncol=ncol(plvalues),dimnames=dimnames(plvalues))
  # qvalue correction
  for (i in 1:ncol(plvalues)) {
    tqs <- tryCatch(qvalue(na.omit(plvalues[,i]))$qvalues, 
                    error = function(e) NULL)
    #print(tqs)
    if (length(tqs) >0) {
      qvalues[names(tqs),i] <- tqs
    }
    else {
      qvalues[names(tqs),i] <- NA
    }
    
  }
  return(list(plvalues=qvalues,Sds=sqrt(lm.bayes$s2.post)))
}

#' Calculate "biological homogeneity index" from GO terms and uniprot accessn names 
#' @export
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

#' Wrapper for statistical testing
#' @importFrom matrixStats rowSds
#' @importFrom shiny validate
#' @export
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

#' Wrapper for estimation of cluster number
#' @import limma
#' @importFrom shiny getDefaultReactiveDomain incProgress
#' @importFrom matrixStats rowMaxs
#' @export
estimClustNum<- function(dat, maxClust=25, cores=1) {
  ClustInd<-matrix(NA,nrow=maxClust,ncol=6)
  tData <- dat[,1:(ncol(dat)-1)]
  colnames(tData)<-NULL
  # We do not filter anymore for NAs, and we do not impute
  # PExpr <- new("ExpressionSet",expr=as.matrix(tData))
  # PExpr.r <- filter.NA(PExpr, thres = 0.25)
  # PExpr <- fill.NA(PExpr.r,mode = "mean")
  # tmp <- filter.std(PExpr,min.std=0,visu=F)
  
  # Standardise
  # PExpr2 <- standardise(PExpr)
  tData <- t(scale(t(tData)))
  sds <- dat[rownames(tData), ncol(dat)]
  
  multiOut <- lapply(3:maxClust,function(x) {    
    if (!is.null(getDefaultReactiveDomain())) {
      incProgress(1, detail = paste("Running cluster number",x))
    } else {
      print(paste("Running cluster number",x))
    }
    clustout <- ClustComp(tData,NClust=x,Sds=sds,NSs=16, cores)
    c(clustout$indices,sum(rowMaxs(clustout$Bestcl$membership)>0.5),
      sum(rowMaxs(clustout$Bestcl2$membership)>0.5))
  })
  for (NClust in 3:maxClust) 
    ClustInd[NClust,] <- multiOut[[NClust-2]]
  
  dmindist <- c(which.max(ClustInd[3:(maxClust-2),1]-ClustInd[4:(maxClust-1),1])+2,
                which.max(ClustInd[3:(maxClust-2),3]-ClustInd[4:(maxClust-1),3])+2)
  dxiebeni <- c(which.min(ClustInd[3:(maxClust-1),2])+2,
                which.min(ClustInd[3:(maxClust-1),4])+2)  
  print(dmindist)
  # graphics.off()
  #  plot.new() ## clean up device
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

#' Wrapper for clustering
#' @importFrom shiny getDefaultReactiveDomain incProgress
#' @importFrom matrixStats rowMaxs
#' @export
runClustWrapper <- function(dat, NClust, proteins=NULL, VSClust=T, cores) {
  tData <- dat[,1:(ncol(dat)-1)]
  
  # We do not filter and impute anymore  
  # PExpr.r <- filter.NA(PExpr, thres = 0.25)
  # PExpr <- fill.NA(PExpr.r,mode = "mean")
  # tmp <- filter.std(PExpr,min.std=0,visu=F)
  #PExpr <- standardise(PExpr)
  
  #Standardize
  tData <- t(scale(t(tData)))
  
  
  clustout <- ClustComp(tData,NClust=NClust,Sds=dat[rownames(tData),ncol(dat)],NSs=16, cores)
  if (VSClust) {
    Bestcl <- clustout$Bestcl
  } else {
    Bestcl <- clustout$Bestcl2
  }
  Bestcl <- SwitchOrder(Bestcl,NClust)
  
  # sorting for membership values (globally)
  Bestcl$cluster <- Bestcl$cluster[order(rowMaxs(Bestcl$membership,na.rm=T))]
  Bestcl$membership <- Bestcl$membership[order(rowMaxs(Bestcl$membership,na.rm=T)),]
  tData <- tData[names(Bestcl$cluster),]
  
  if (!is.null(getDefaultReactiveDomain()))
    incProgress(0.7, detail = paste("Plotting",NClust))
  
  # graphics.off() ## clean up device
  par(lwd=0.25)
  oldmar <- par("mar")
  par(mar=c(2,2,3,3),mgp=c(2,1,0))
  par(mar=par("mar")/max(1,NClust/20))
  
  ## Plot results
  vs_filename <- NULL
  if (VSClust) {
    vs_filename <- "VSClust_New_Clusters.pdf"
  } else {
    vs_filename <- "VSClust_Std_Clusters.pdf"
  }
  mfuzz.plot(tData,cl=Bestcl,mfrow=c(round(sqrt(NClust)),ceiling(sqrt(NClust))),min.mem=0.5,colo="fancy", filename = vs_filename)
  mfuzz.plot(tData,cl=Bestcl,mfrow=c(round(sqrt(NClust)),ceiling(sqrt(NClust))),min.mem=0.5,colo="fancy")
  p <- recordPlot()
  # par(lwd=1,mar=oldmar)
  
  colnames(Bestcl$membership) <- paste("membership of cluster",colnames(Bestcl$membership))
  outFileClust <- tData
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
  Out <- list(dat=tData, Bestcl=Bestcl, p=p, outFileClust=outFileClust, ClustInd=ClustInd)
  return(Out)
}

# Wrapper for functional enrichment
#' @importFrom clusterProfiler compareCluster
#' @importFrom matrixStats rowMaxs
#' @export
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
  if (!is.null(getDefaultReactiveDomain()))
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

#' Plotting clustering results into multiple figure panels, adopted from the MFuzz package
#' @param filename for writing into pdf. Will write on screen when using NA
#' @export

mfuzz.plot <- function (dat, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels, 
                        filename=NA,xlab="Time",ylab="Expression changes") 
{
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1
  colorindex <- integer(dim(dat)[[1]])
  if (missing(colo)) {
    colo <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", 
              "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", 
              "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", 
              "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", 
              "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", 
              "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", 
              "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", 
              "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", 
              "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", 
              "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", 
              "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", 
              "#FF0060", "#FF0048", "#FF0030", "#FF0018")
  }    else {
    if (colo == "fancy") {
      fancy.blue <- c(c(255:0), rep(0, length(c(255:0))), 
                      rep(0, length(c(255:150))))
      fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
      fancy.red <- c(c(0:255), rep(255, length(c(255:0))), 
                     c(255:150))
      colo <- rgb(b = fancy.blue/255, g = fancy.green/255, 
                  r = fancy.red/255)
    }
  }
  colorseq <- seq(0, 1, length = length(colo))
  for (j in 1:max(clusterindex)) {
    tmp <- dat[clusterindex == j, ]
    tmpmem <- memship[clusterindex == j, j]
    if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0) {
      if (!is.na(filename)) {
        pdf(filename, height=3*mfrow[1],width=3*mfrow[2])
      }
      par(mfrow = mfrow, cex=0.5)
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      plot.default(x = NA, xlim = c(1, dim(dat)[[2]]), 
                   ylim = c(ymin, ymax), xlab = xlab, ylab = ylab, 
                   main = paste("Cluster", j), axes = FALSE)
      if (missing(time.labels)) {
        axis(1, 1:dim(dat)[[2]], c(1:dim(dat)[[2]]))
        axis(2)
      }
      else {
        axis(1, 1:dim(dat)[[2]], time.labels)
        axis(2)
      }
    }
    else {
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      plot.default(x = NA, xlim = c(1, dim(dat)[[2]]), 
                   ylim = c(ymin, ymax), xlab = xlab, ylab = ylab, 
                   main = paste("Cluster", j), axes = FALSE)
      if (missing(time.labels)) {
        axis(1, 1:dim(dat)[[2]], c(1:dim(dat)[[2]]))
        axis(2)
      }
      else {
        axis(1, 1:dim(dat)[[2]], time.labels)
        axis(2)
      }
    }
    if (!(sum(clusterindex == j) == 0)) {
      for (jj in 1:(length(colorseq) - 1)) {
        tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= 
                     colorseq[jj + 1])

        if (sum(tmpcol, na.rm=T) > 0) {
          tmpind <- which(tmpcol)
          for (k in 1:length(tmpind)) {
            lines(tmp[tmpind[k], ], col = colo[jj])
          }
        }
      }
    }
  }
  if (!is.na(filename))   dev.off()
}


