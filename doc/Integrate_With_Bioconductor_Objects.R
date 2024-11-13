## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----env, message = FALSE, warning = FALSE, echo = FALSE----------------------
library("vsclust")
library("MultiAssayExperiment")

## ----eval=FALSE---------------------------------------------------------------
# # uncomment in case you have not installed vsclust yet
# #if (!require("BiocManager", quietly = TRUE))
# #    install.packages("BiocManager")
# #BiocManager::install("vsclust")

## -----------------------------------------------------------------------------
#### Input parameters, only read when now parameter file was provided #####
## All principal parameters for running VSClust can be defined as in the 
## shiny app at computproteomics.bmb.sdu.dk/Apps/VSClust
# name of study
Experiment <- "miniACC" 
# Paired or unpaired statistical tests when carrying out LIMMA for 
# statistical testing
isPaired <- FALSE
# Number of threads to accelerate the calculation (use 1 in doubt)
cores <- 1 

# If 0 (default), then automatically estimate the cluster number for the 
# vsclust run from the Minimum Centroid Distance
PreSetNumClustVSClust <- 0 
# If 0 (default), then automatically estimate the cluster number for the 
# original fuzzy c-means from the Minimum Centroid Distance
PreSetNumClustStand <- 0 

# max. number of clusters when estimating the number of clusters. 
# Higher numbers can drastically extend the computation time.
maxClust <- 10 

## ----fig.width = 12-----------------------------------------------------------
data(miniACC, package="MultiAssayExperiment")
# log-transformation and remove of -Inf values
logminiACC <- log2(assays(miniACC)$RNASeq2GeneNorm)
logminiACC[!is.finite(logminiACC)] <- NA
# normalize to median
logminiACC <- t(t(logminiACC) - apply(logminiACC, 2, median, na.rm=TRUE))

miniACC2 <- c(miniACC, log2rnaseq = logminiACC, mapFrom=1L)

boxplot(logminiACC)


#### running statistical analysis and estimation of individual variances
statOut <- PrepareSEForVSClust(miniACC2, "log2rnaseq", 
                               coldatname = "OncoSign", 
                               isPaired=isPaired, isStat=TRUE)


## ----fig.width = 12-----------------------------------------------------------

#### Estimate number of clusters with maxClust as maximum number clusters to run 
#### the estimation with
ClustInd <- estimClustNum(statOut$dat, maxClust=maxClust, cores=cores)

#### Use estimate cluster number or use own
if (PreSetNumClustVSClust == 0)
  PreSetNumClustVSClust <- optimalClustNum(ClustInd)
if (PreSetNumClustStand == 0)
  PreSetNumClustStand <- optimalClustNum(ClustInd, method="FCM")
#### Visualize
  estimClust.plot(ClustInd)




## -----------------------------------------------------------------------------


#### Run clustering (VSClust and standard fcm clustering
ClustOut <- runClustWrapper(statOut$dat, 
                            PreSetNumClustVSClust, NULL, 
                            VSClust=TRUE,
                            cores=cores)
Bestcl <- ClustOut$Bestcl
VSClust_cl <- Bestcl


## -----------------------------------------------------------------------------
sessionInfo()

