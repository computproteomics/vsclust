########### functions for fuzzyVarM

cores <- 4


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


ClustComp <- function(tData,NSs=10,NClust=NClust,Sds=Sds) {
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
  
  cls <- mclapply(1:NSs, function(x) cmeans(exprs(PExpr2),NClust,m=m,verbose=F,iter.max=1000), mc.cores=cores)
  Bestcl <- cls[[which.min(lapply(cls,function(x) x$withinerror))]]
  # for (cc in 1:length(cls)) {
  #   print(paste("New error:",cls[[cc]]$withinerror, "iter:",cls[[cc]]$iter))
  # }
  cls <- mclapply(1:NSs, function(x) cmeans(exprs(PExpr2),NClust,m=mm,verbose=F,iter.max=1000), mc.cores=cores)
  Bestcl2 <- cls[[which.min(lapply(cls,function(x) x$withinerror))]]
  # for (cc in 1:length(cls)) {
  #   print(paste("Std error:",cls[[cc]]$withinerror, "iter:",cls[[cc]]$iter))
  # }
  # 
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
    tqs <- qvalue(na.omit(plvalues[,i]))$qvalues
    qvalues[names(tqs),i] <- tqs
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
    tqs <- qvalue(na.omit(plvalues[,i]))$qvalues
    qvalues[names(tqs),i] <- tqs
  }
  return(list(plvalues=qvalues,Sds=sqrt(lm.bayes$s2.post)))
}



## overwrite function to set timeout limit higher
enrichDAVID2 <- function (gene, idType = "ENTREZ_GENE_ID", listType = "Gene", 
                          minGSSize = 5, annotation = "GOTERM_BP_ALL", pvalueCutoff = 0.05, 
                          pAdjustMethod = "BH", qvalueCutoff = 0.2, species = NA, david.user = "veits@bmb.sdu.dk") 
{
  Count <- List.Total <- Pop.Hits <- Pop.Total <- NULL
  pAdjustMethod <- match.arg(pAdjustMethod, c("bonferroni", 
                                              "BH"))
  idType <- match.arg(idType, c("AFFYMETRIX_3PRIME_IVT_ID",
                                "AFFYMETRIX_EXON_GENE_ID", "AFFYMETRIX_SNP_ID", "AGILENT_CHIP_ID",
                                "AGILENT_ID", "AGILENT_OLIGO_ID", "ENSEMBL_GENE_ID",
                                "ENSEMBL_TRANSCRIPT_ID", "ENTREZ_GENE_ID", "GENOMIC_GI_ACCESSION",
                                "GENPEPT_ACCESSION", "ILLUMINA_ID", "IPI_ID", "MGI_ID",
                                "OFFICIAL_GENE_SYMBOL", "PFAM_ID", "PIR_ID", "PROTEIN_GI_ACCESSION",
                                "REFSEQ_GENOMIC", "REFSEQ_MRNA", "REFSEQ_PROTEIN", "REFSEQ_RNA",
                                "RGD_ID", "SGD_ID", "TAIR_ID", "UCSC_GENE_ID", "UNIGENE",
                                "UNIPROT_ACCESSION", "UNIPROT_ID", "UNIREF100_ID", "WORMBASE_GENE_ID",
                                "WORMPEP_ID", "ZFIN_ID"))
  david <- DAVIDWebService$new(email = david.user,url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  setTimeOut(david,1000000)
  david.res <- addList(david, gene, idType = idType, listName = "clusterProfiler", 
                       listType = listType)
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
  ID <- term.df[, 1]
  Description <- term.df[, 2]
  GeneRatio <- with(x, paste(Count, List.Total, sep = "/"))
  BgRatio <- with(x, paste(Pop.Hits, Pop.Total, sep = "/"))
  Over <- data.frame(ID = ID, Description = Description, GeneRatio = GeneRatio, 
                     BgRatio = BgRatio, pvalue = x$PValue)
  row.names(Over) <- ID
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
  gc <- strsplit(Over$geneID, "/")
  names(gc) <- Over$ID
  new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff, 
      pAdjustMethod = pAdjustMethod, organism = org, ontology = as.character(x$Category[1]), 
      gene = as.character(gene), geneInCategory = gc)
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
