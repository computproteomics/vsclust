library(matrixStats)
library(Mfuzz)
library(limma)
library(parallel)
library(qvalue)
library(e1071FuzzVec)
library(shinyjs)
library(clusterProfiler)
#library(RDAVIDWebService)
source("FcmClustPEst.R")
source("mfuzz.plotpdf.R")
source("HelperFuncs.R")

dat <- read.csv("../../AllProteins.csv",row.names=1)
dat <- dat[1:nrow(dat),]
Accs <- list()
for (c in 1:5) {
  Accs[[c]] <- sample(rownames(dat),100)
  Accs[[c]] <- Accs[[c]][Accs[[c]]!=""]
  if (length(Accs[[c]])>0) {
    if (length(Accs[[c]])>1) {
      tdat <- (dat[Accs[[c]],])
    } else {
      tdat <- t(dat[Accs[[c]],])
    }
    Accs[[c]] <- sub("-[0-9]","",Accs[[c]])
  } else {
    tdat <- t(rep(NA,ncol(dat)+1))
  }
}
names(Accs) <- paste("Cluster",1:length(Accs))
Accs <- lapply(Accs,function(x) unique(ifelse(is.na(x),"B3",x)))
Accs <- Accs[lapply(Accs,length)>0]
x <- compareCluster(Accs, fun="enrichDAVID2", annotation="GOTERM_BP_ALL",
                    idType="UNIPROT_ACCESSION",
                    listType="Gene", david.user = "veits@bmb.sdu.dk")
x@compareClusterResult <- cbind(x@compareClusterResult,log10padval=log10(x@compareClusterResult$p.adjust))
BHI <- calcBHI(Accs,x)
y <- new("compareClusterResult",compareClusterResult=x@compareClusterResult)
if (length(unique(y@compareClusterResult$ID)) > 20) {
  print("Reducing number of DAVID results")
  y@compareClusterResult <- y@compareClusterResult[
    order(y@compareClusterResult$p.adjust)[1:20],]
  # print(x@compareClusterResult)
}
plot(y,title=paste("BHI:",BHI),showCategory=1000,colorBy="log10padval",font.size=10)
