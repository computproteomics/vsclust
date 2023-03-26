#!/usr/bin/env Rscript

## Initialization ### 
library(vsclust)
library(yaml)
library(shiny)
library(clusterProfiler)
library(matrixStats)

############ Command-line wrapper for VSClust

######## reading parameter file
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  stop("VSClust can be run either by defining the parameters in the script or by specifying a parameter when the calling runVSClust.R.
    For setting the parameters, see vsclust.yml in the main folder as example and descriptions of the parameter values.", call.=FALSE)
} else if (length(args) == 1) {
  parfile <- args[1]
  print(paste("Reading",parfile))
  
  ## read all parameters from yaml file + checks
  pars <- yaml.load_file(parfile)
  Experiment <- pars$Experiment
  
  NumReps <- pars$NumReps
  if (!is.integer(NumReps))
    stop("Number of replicates is not a positive integer number", call.=FALSE)
  if (NumReps < 0)
    stop("Number of replicates is not a positive integer number", call.=FALSE)
  NumCond <- pars$NumCond
  if (!is.integer(NumCond)) 
    stop("Number of conditions is not a positive integer number", call.=FALSE)
  if (NumCond < 0) 
    stop("Number of conditions is not a positive integer number", call.=FALSE)
  isPaired <- pars$paired
  if (!is.logical(isPaired)) 
    stop("paired parameter should be 'true' or 'false'")
  isStat <- pars$stat
  if (!is.logical(isStat)) 
    stop("stat parameter should be 'true' or 'false'")
  infile <- pars$infile
  if (!file.exists(infile))
    stop(paste0("input file ",infile," does not exist"), call.=FALSE)
  protnames <- pars$secondcol
  if (!is.logical(protnames)) 
    stop("secondcol parameter should be 'true' or 'false'")
  is_header <- pars$is_header
  if (!is.logical(protnames)) 
    stop("is_header parameter should be 'true' or 'false'")
  cores <- pars$cores
  if (!is.integer(cores)) 
    stop("Parameter cores is not a positive integer number", call.=FALSE)
  if (cores < 1 ) 
    stop("Parameter cores is not a positive integer number > 0", call.=FALSE)
  PreSetNumClustVSClust <- pars$PreSetNumClustVSClust
  if (!is.integer(PreSetNumClustVSClust)) 
    stop("Parameter PreSetNumClustVSClust is not a positive integer number", call.=FALSE)
  if (PreSetNumClustVSClust < 1 & PreSetNumClustVSClust != 0) 
    stop("Parameter PreSetNumClustVSClust is not a positive integer number or 0", call.=FALSE)
  PreSetNumClustStand <- pars$PreSetNumClustStand
  if (!is.integer(PreSetNumClustStand)) 
    stop("Parameter PreSetNumClustStand is not a positive integer number", call.=FALSE)
  if (PreSetNumClustStand < 1 & PreSetNumClustStand != 0) 
    stop("Parameter PreSetNumClustStand is not a positive integer number or 0", call.=FALSE)
  maxClust <- pars$maxClust
  if (!is.integer(maxClust)) 
    stop("Parameter maxClust is not a positive integer number", call.=FALSE)
  if (maxClust < 3) 
    stop("Parameter maxClust should be at least 3", call.=FALSE)
  
} else {
  
  
  #### Input parameters, only read when now parameter file was provided #####
  # All principal parameters for running VSClust can be defined as in the shiny app at computproteomics.bmb.sdu.dk/Apps/VSClust
  Experiment <- "ProtExample" ## name of study
  NumReps <- 3###886 ## Number of replicates per
  NumCond <- 4###12 ## Number of different experimental conditions (e.g. time points)
  isPaired <- F ## Paired or unpaired statistical tests
  isStat <- T ## Set to F when no replicates but last column with individual standard deviations
  infile <- system.file('extdat/ProtExample.csv', package='vsclust')##"/path/MyData.csv" ## Input filename
  protnames <- F ## Low-level data (e.g. probes of transcripts or peptides)
  is_header <- T ## File contains one-line header
  cores <- 4#4 # Number of cores to use ## 1 is for windows
  
  PreSetNumClustVSClust <- 0 # If 0, then automatically take the one from Minimum Centroid Distance
  PreSetNumClustStand <- 0 # If 0, then automatically take the one from Minimum Centroid Distance
  
  
  maxClust <- 20 ## max. number of clusters when estimating the number of clusters
}


## reading helper functions
# need to change to the source path and back
currPath <- getwd()
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
cat(paste0("Getting R functions from files located in ", script.basename),"\n")
setwd(script.basename)
setwd(currPath)


#### File readin

dat <- read.csv(infile,row.names=1,header=is_header)    
dat <- dat[rownames(dat)!="",]
proteins <- NULL
if(protnames) {
  proteins <- dat[,1]
  dat <- dat[,2:ncol(dat)]
  names(proteins) <- rownames(dat)
}

#### In case you need grouping is on putting replicates in adjacent columns -> reorganize order
dat <- dat[,rep(0:(NumCond-1),NumReps)*NumReps+rep(1:(NumReps), each=NumCond)]



#### running statistical analysis and estimation of individual variances
statOut <- PrepareForVSClust(dat, NumReps, NumCond, isPaired, isStat)

dat <- statOut$dat
Sds <- dat[,ncol(dat)]
print(paste("Features:",nrow(dat),"<br/>Missing values:",
            sum(is.na(dat)),"<br/>Median standard deviations:",
            round(median(Sds,na.rm=T),digits=3)))

## Create pdf-figure of validation indices "minimum centroid distance" and "Xie-Beni index"
write.csv(statOut$statFileOut,paste("",Experiment,"statFileOut.csv",sep=""))

#### Estimate number of clusters with maxClust as maximum number clusters to test for
## Write output into file 
pdf(paste(Experiment,"EstimatedClustNumber.pdf", sep=""),height=6,width=15)
ClustInd <- estimClustNum(dat, maxClust, cores=cores)
estimClust.plot(ClustInd)
dev.off()

#### Use estimate cluster number or use own
if (PreSetNumClustVSClust == 0)
  PreSetNumClustVSClust <- optimalClustNum(ClustInd)
if (PreSetNumClustStand == 0)
  PreSetNumClustStand <- optimalClustNum(ClustInd, method="FCM")



#### Run clustering (VSClust and standard fcm clustering
## Write pdf-figure of clusters
pdf(paste(Experiment,"FCMVarMResults", Sys.Date(), ".pdf", sep=""))
ClustOut <- runClustWrapper(dat, PreSetNumClustVSClust, proteins, VSClust=T, cores=cores)
dev.off() 

Bestcl <- ClustOut$Bestcl
ClustOut$p
## Write clustering results (VSClust)
write.csv(data.frame(cluster=Bestcl$cluster,ClustOut$outFileClust,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                     Bestcl$membership), paste(Experiment, "FCMVarMResults", Sys.Date(), ".csv", sep=""))
## Write coordinates of cluster centroids
write.csv(Bestcl$centers, paste(Experiment,"FCMVarMResultsCentroids", Sys.Date(), ".csv", sep=""))

## Write pdf-figure of clusters
pdf(paste(Experiment,"FCMStdResults", Sys.Date(), ".pdf", sep=""))
ClustOut <- runClustWrapper(dat, PreSetNumClustStand, proteins, VSClust=F, cores=cores)
dev.off()

Bestcl <- ClustOut$Bestcl
ClustOut$p
## Write clustering results (standard fcm)
write.csv(data.frame(cluster=Bestcl$cluster,ClustOut$outFileClust,isClusterMember=rowMaxs(Bestcl$membership)>0.5,maxMembership=rowMaxs(Bestcl$membership),
                     Bestcl$membership), paste(Experiment, "FCMResults", Sys.Date(), ".csv", sep=""))
## Write coordinates of cluster centroids
write.csv(Bestcl$centers, paste(Experiment,"FCMResultsCentroids", Sys.Date(), ".csv", sep=""))


