## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----env, message = FALSE, warning = FALSE, echo = FALSE----------------------
require(clusterProfiler)
require(matrixStats)

## ----eval=FALSE---------------------------------------------------------------
# # Uncomment in case you have not installed vsclust yet
# #if (!require("BiocManager", quietly = TRUE))
# #    install.packages("BiocManager")
# #BiocManager::install("vsclust")
# library(vsclust)

## -----------------------------------------------------------------------------
#### Input parameters, only read when now parameter file was provided
## All principal parameters for running VSClust can be defined as in the 
## shinyapp at computproteomics.bmb.sdu.dk/Apps/VSClust 
# name of study
Experiment <- "ProtExample" 
# Number of replicates/sample per different experimental condition (sample 
# type)
NumReps <- 3  
# Number of different experimental conditions (e.g. time points or sample 
# types)
NumCond <- 4  
# Paired or unpaired statistical tests when carrying out LIMMA for 
# statistical testing
isPaired <- FALSE
# Number of threads to accelerate the calculation (use 1 in doubt)
cores <- 1 

# If 0 (default), then automatically estimate the cluster number for the 
# vsclust 
# run from the Minimum Centroid Distance
PreSetNumClustVSClust <- 0 
# If 0 (default), then automatically estimate the cluster number for the 
# original fuzzy c-means from the Minimum Centroid Distance
PreSetNumClustStand <- 0 

# max. number of clusters when estimating the number of clusters. Higher 
# numbers can drastically extend the computation time.
maxClust <- 10 

## ----fig.width = 12-----------------------------------------------------------

data(protein_expressions)
dat <- protein_expressions

#### running statistical analysis and estimation of individual variances
statOut <- PrepareForVSClust(dat, NumReps, NumCond, isPaired, TRUE)

dat <- statOut$dat
Sds <- dat[,ncol(dat)]
cat(paste("Features:",nrow(dat),"\nMissing values:",
            sum(is.na(dat)),"\nMedian standard deviations:",
            round(median(Sds,na.rm=TRUE),digits=3)))

## Write output into file 
write.csv(statOut$statFileOut,
          paste("",Experiment,"statFileOut.csv",sep=""))


## ----fig.width = 12-----------------------------------------------------------

#### Estimate number of clusters with maxClust as maximum number clusters 
#### to run the estimation with
ClustInd <- estimClustNum(dat, maxClust=maxClust, scaling="standardize", cores=cores)

#### Use estimate cluster number or use own
if (PreSetNumClustVSClust == 0)
  PreSetNumClustVSClust <- optimalClustNum(ClustInd)
if (PreSetNumClustStand == 0)
  PreSetNumClustStand <- optimalClustNum(ClustInd, method="FCM")
#### Visualize
  estimClust.plot(ClustInd)


## -----------------------------------------------------------------------------
#### Run clustering (VSClust and standard fcm clustering
ClustOut <- runClustWrapper(dat, 
                            PreSetNumClustVSClust, 
                            NULL, 
                            VSClust=TRUE, 
                            scaling="standardize",
                            cores=cores)
Bestcl <- ClustOut$Bestcl
VSClust_cl <- Bestcl
#ClustOut$p
## Write clustering results (VSClust)
write.csv(data.frame(cluster=Bestcl$cluster,
                     ClustOut$outFileClust,
                     isClusterMember=rowMaxs(Bestcl$membership)>0.5,
                     maxMembership=rowMaxs(Bestcl$membership),
                     Bestcl$membership), 
          paste(Experiment, 
                "FCMVarMResults", 
                Sys.Date(), 
                ".csv", 
                sep=""))
## Write coordinates of cluster centroids
write.csv(Bestcl$centers, 
          paste(Experiment,
                "FCMVarMResultsCentroids",
                Sys.Date(), 
                ".csv", 
                sep=""))

## -----------------------------------------------------------------------------
ClustOut <- runClustWrapper(dat, PreSetNumClustStand, NULL, VSClust=FALSE, 
                            scaling="standardize", cores=cores)
Bestcl <- ClustOut$Bestcl
## Write clustering results (standard fcm)
write.csv(data.frame(cluster=Bestcl$cluster,
                     ClustOut$outFileClust,
                     isClusterMember=rowMaxs(Bestcl$membership)>0.5,
                     maxMembership=rowMaxs(Bestcl$membership),
                     Bestcl$membership), 
          paste(Experiment, 
                "FCMResults", 
                Sys.Date(), 
                ".csv", 
                sep=""))
## Write coordinates of cluster centroids
write.csv(Bestcl$centers, paste(Experiment,
                                "FCMResultsCentroids", 
                                Sys.Date(),
                                ".csv", 
                                sep=""))

## -----------------------------------------------------------------------------
sessionInfo()

