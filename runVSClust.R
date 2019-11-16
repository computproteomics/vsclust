############ Command-line wrapper for VSClust

#### Input parameters #####
# All principal parameters for running VSClust can be defined as in the shiny app at computproteomics.bmb.sdu.dk/Apps/VSClust

Experiment <- "ProtExample" ## name of study
NumReps <- 3###886 ## Number of replicates per
NumCond <- 4###12 ## Number of different experimental conditions (e.g. time points)
isPaired <- F ## Paired or unpaired statistical tests
isStat <- T ## Set to F when no replicates but additional columns with individual standard deviations
infile <- 'ProtExample.csv'##"/path/MyData.csv" ## Input filename
protnames <- F ## Low-level data (e.g. probes of transcripts or peptides)
is_header <- T ## File contains one-line header
cores <- 4#4 # Number of cores to use ## 1 is for windows

PreSetNumClustVSClust <- 0 # If 0, then automatically take the one from Minimum Centroid Distance
PreSetNumClustStand <- 0 # If 0, then automatically take the one from Minimum Centroid Distance


maxClust <- 20 ## max. number of clusters when estimating the number of clusters


## packages
library(matrixStats)
library(Mfuzz)
library(limma)
library(qvalue)
## HAS TO BE THE MODIFIED ONE!
require(e1071FuzzVec)
library(shiny)
library(clusterProfiler)
library(RDAVIDWebService)
source("FcmClustPEst.R")
source("mfuzz.plotpdf.R")
source("HelperFuncs.R")


#### File readin

dat <- read.csv(infile,row.names=1,header=is_header)    
dat <- dat[rownames(dat)!="",]
proteins <- NULL
if(protnames) {
  proteins <- dat[,1]
  dat <- dat[,2:ncol(dat)]
  names(proteins) <- rownames(dat)
}


#### running statistical analysis and estimation of individual variances
statOut <- statWrapper(dat, NumReps, NumCond, isPaired, isStat)

dat <- statOut$dat
Sds <- dat[,ncol(dat)]
print(paste("Features:",nrow(dat),"<br/>Missing values:",
            sum(is.na(dat)),"<br/>Median standard deviations:",
            round(median(Sds,na.rm=T),digits=3)))

## Write output into file 
write.csv(statOut$statFileOut,paste("",Experiment,"statFileOut.csv",sep=""))

#### Estimate number of clusters with maxClust as maximum number clusters to test for
clustNumOut <- estimClustNum(dat, maxClust, cores)

#### Use estimate cluster number or use own
if (PreSetNumClustVSClust == 0)
  PreSetNumClustVSClust <- clustNumOut$numclust
if (PreSetNumClustStand == 0)
  PreSetNumClustStand <- clustNumOut$numclust

## Create pdf-figure of validation indices "minimum centroid distance" and "Xie-Beni index"
pdf(paste(Experiment,"EstimatedClustNumber.pdf", sep=""),height=6,width=15)
print(clustNumOut$p)
dev.off()  


#### Run clustering (VSClust and standard fcm clustering
ClustOut <- runClustWrapper(dat, PreSetNumClustVSClust, proteins, VSClust=T, cores)
Bestcl <- ClustOut$Bestcl
ClustOut$p
## Write clustering results (VSClust)
write.csv(data.frame(cluster=Bestcl$cluster,ClustOut$outFileClust,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                     Bestcl$membership), paste(Experiment, "FCMVarMResults", Sys.Date(), ".csv", sep=""))
## Write coordinates of cluster centroids
write.csv(Bestcl$centers, paste(Experiment,"FCMVarMResultsCentroids", Sys.Date(), ".csv", sep=""))

## Write pdf-figure of clusters
pdf(paste(Experiment,"FCMVarMResults", Sys.Date(), ".pdf", sep=""),height=5*round(sqrt(PreSetNumClustVSClust)),width=5*ceiling(sqrt(PreSetNumClustVSClust)))
print(ClustOut$p)
dev.off()
print(ClustOut$ClustInd)

ClustOut <- runClustWrapper(dat, PreSetNumClustStand, proteins, VSClust=F, cores)
Bestcl <- ClustOut$Bestcl
ClustOut$p
## Write clustering results (standard fcm)
write.csv(data.frame(cluster=Bestcl$cluster,ClustOut$outFileClust,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                     Bestcl$membership), paste(Experiment, "FCMResults", Sys.Date(), ".csv", sep=""))
## Write coordinates of cluster centroids
write.csv(Bestcl$centers, paste(Experiment,"FCMResultsCentroids", Sys.Date(), ".csv", sep=""))

## Write pdf-figure of clusters
pdf(paste(Experiment,"FCMResults", Sys.Date(), ".pdf", sep=""),height=5*round(sqrt(PreSetNumClustStand)),width=5*ceiling(sqrt(PreSetNumClustStand)))
print(ClustOut$p)
dev.off()
print(ClustOut$ClustInd)





