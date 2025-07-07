#' Functions for running VSClust analysis
#'
#' Wrapper for statistical analysis
#'
#' Prepare data for running vsclust clustering.
#' This includes visualization running the functions for the principal component
#' analysis and its visualization, statistical testing with LIMMA, as well as
#' scaling and filtering of missing values
#' @param dat matrix or data frame of numerical data. Columns are samples.
#' Replicates are grouped (i.e. A1, B1, C1, A2, B2, C2) when letters denote
#' conditions and numbers the replicates. In case of `isStat=FALSE`, you need a
#' last column for the standard deviations
#' @param NumReps Number replicates in the data
#' @param NumCond Number of different experimental conditions. The total number
#' of columns needs to be NumReps*NumCond
#' @param isPaired Boolean for running paired or unpaired statistical tests
#' @param isStat Boolean for whether to run statistical test or each column
#' corresponds to a different experimental conditions. Then this function reads
#' feature standard deviations from data frame from the last column
#' @return list with the items `dat` (data matrix of features averaged over
#' replicates and last column with their standard deviations), `qvals` FDRs from
#' the statistical tests (each conditions versus the first), `StatFileOut` all
#' of before for saving in file
#' @examples
#' data <- matrix(rnorm(2000), nrow=200)
#' stats <- PrepareForVSClust(data, 5, 2, isStat=TRUE)
#'
#' @import stats
#' @importFrom matrixStats rowSds
#' @importFrom shiny validate
#' @export
#' @references
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering
#' of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi:
#' 10.1093/bioinformatics/bty224. PMID: 29635359.
#'
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and
#' the Quantitative Analysis of Protein Complexes. Methods Mol Biol.
#' 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#'
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters
#' for fuzzy c-means cluster analysis. Bioinformatics. 2010
#' Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29.
#' PMID: 20880957.
PrepareForVSClust <-
    function(dat, NumReps, NumCond, isPaired = FALSE, isStat) {
        qvals <- statFileOut <- Sds <- NULL
        tdat <- NULL
        
        # convert to matrix
        dat <- as.matrix(dat)
        
        # Run statistical testing
        if (isStat) {
            if (ncol(dat) != NumReps * NumCond)
                stop("Number of data columns must correspond to product of conditions
             and replicates!")
            if (isPaired) {
                ttt <- SignAnalysisPaired(dat, NumCond, NumReps)
            } else {
                ttt <- SignAnalysis(dat, NumCond, NumReps)
            }
            
            Sds <- ttt$Sds
            qvals <- ttt$qvalues
            colnames(qvals) <-
                paste("qvalue ", LETTERS702[2:(NumCond)], "vsA", sep = "")
            
            tdat <- averageCond(dat, NumReps, NumCond)
            
        } else {
            Sds <- dat[, ncol(dat)]
            tdat <- dat[, seq_len(ncol(dat) - 1)]
            NumReps <- 1
            NumCond <- ncol(dat) - 1
            dat <- tdat
        }
        
        if (isStat) {
            statFileOut <- cbind(tdat, Sds, qvals)
        } else {
            statFileOut <- cbind(tdat, Sds)
        }
        
        pcaWithVar(dat, NumReps, NumCond, Sds / rowSds(tdat, na.rm = TRUE))
        
        ## Preparing output
        Out <-
            list(dat = cbind(tdat, Sds),
                 qvals = qvals,
                 statFileOut = statFileOut)
        Out
        
    }

#' Wrapper for statistical analysis for SummarizedExperiment object
#'
#' Prepare data for running vsclust clustering.
#' This includes visualization running the functions for the principal component
#' analysis and its visualization, statistical testing with LIMMA, as well as
#' scaling and filtering of missing values
#' @param se SummarizedExperiment object
#' @param assayname Sample in SummarizedExperiment object
#' @param coldatname Column in colData for extracting replicates
#' @param isPaired Boolean for running paired or unpaired statistical tests
#' @param isStat Boolean for whether to run statistical test or each column
#' corresponds to a different experimental conditions. Then this function reads
#' feature standard deviations from data frame from the last column
#' @return list with the items `dat` (data matrix of features averaged over
#' replicates and last column with their standard deviations), `qvals` FDRs from
#' the statistical tests (each conditions versus the first), `StatFileOut` all
#' of before for saving in file, `NumReps` number of replicates and `NumCond` 
#' number of different experimental conditions
#' @examples
#' data(miniACC, package="MultiAssayExperiment")
#' 
#' stats <- PrepareSEForVSClust(miniACC, coldatname="COC", isStat=TRUE)
#'
#' @import stats
#' @importFrom MultiAssayExperiment assay assays sampleMap colData
#' @importFrom matrixStats rowSds
#' @importFrom shiny validate
#' @export
#' @references
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering
#' of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi:
#' 10.1093/bioinformatics/bty224. PMID: 29635359.
#'
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and
#' the Quantitative Analysis of Protein Complexes. Methods Mol Biol.
#' 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#'
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters
#' for fuzzy c-means cluster analysis. Bioinformatics. 2010
#' Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29.
#' PMID: 20880957.
PrepareSEForVSClust <-
    function(se,
             assayname = 1,
             coldatname = NULL,
             isPaired = FALSE,
             isStat) {
        qvals <- statFileOut <- Sds <- NULL
        tdat <- NULL
        
        if (!(class(se) %in% c(
            "SummarizedExperiment",
            "QFeatures",
            "MultiAssayExperiment"
        ))) {
            stop(
                "!! First argument must be a SummarizedExperiment, QFeatures or MultiAssayExperiment objectd"
            )
        }
        
        # convert to matrix
        dat <- assay(se, assayname)
        
        # determine number of conditions and replicates from colData
        NumCond <- NumReps <- 0
        if (!is.null(coldatname)) {
            # change to name if assay is given by index
            if (is.numeric(assayname))
                assayname <- names(se)[assayname]
            sample_names <- sampleMap(se)
            sample_names <- sample_names[sample_names$assay == assayname,]
            rownames(sample_names) <- sample_names$colname
            coldat <- colData(se)[sample_names[colnames(dat), "primary"], coldatname]
            names(coldat) <- colnames(dat)
            NumReps <- max(table(coldat))
            NumCond <- length(unique(coldat))
            message("-- The following categories will be used as experimental 
              conditions:\n",paste(unique(coldat), collapse="\n"))
            if (length(unique(coldat)) < 3)
                stop("!! We need a minimum of three different categories/conditions")
            message("-- Extracted NumReps: ", NumReps, " and NumCond: ", NumCond)
            dat <- balanceData(dat, coldat)
        } else{
            NumCond <- ncol(dat)
            NumReps <- 1
            message("-- No replicates given or no statistical testing, assuming that each
              sample is a different type of sample. Variances will be set to 1")
            dat <- cbind(dat, 1)
        }
        
        # Run statistical testing
        if (isStat) {
            if (ncol(dat) != NumReps * NumCond)
                stop("!! Number of data columns must correspond to product of conditions
             and replicates!")
            if (isPaired) {
                ttt <- SignAnalysisPaired(dat, NumCond, NumReps)
            } else {
                ttt <- SignAnalysis(dat, NumCond, NumReps)
            }
            
            Sds <- ttt$Sds
            qvals <- ttt$qvalues
            colnames(qvals) <-
                paste("qvalue ", LETTERS702[2:(NumCond)], "vsA", sep = "")
            
            tdat <- averageCond(dat, NumReps, NumCond)
            
        } else {
            Sds <- dat[, ncol(dat)]
            tdat <- dat[, seq_len(ncol(dat) - 1)]
            NumReps <- 1
            NumCond <- ncol(dat) - 1
            dat <- tdat
        }
        
        if (isStat) {
            statFileOut <- cbind(tdat, Sds, qvals)
        } else {
            statFileOut <- cbind(tdat, Sds)
        }
        
        pcaWithVar(dat, NumReps, NumCond, Sds / rowSds(tdat, na.rm = TRUE))
        
        ## Preparing output
        Out <-
            list(dat = cbind(tdat, Sds),
                 qvals = qvals,
                 statFileOut = statFileOut, NumReps=NumReps, NumCond=NumCond)
        Out
        
    }


#' Wrapper for estimation of cluster number
#'
#' This runs the clustering for different numbers of clusters, and estimates the
#' most suitable numbers from applying the minimum centroid distance and the Xie
#' Beni index. Multi-threading is used to shorten the computation times.
#' Given the hierarchical structure of many data sets, the resulting numbers are
#' suggestions. Inspection of the here plotted indices help to determine
#' alternative cluster numbers, given by a strong decay of the minimum centroid
#' distance and/or a low value of the Xie Beni index.
#'
#' @param dat matrix of features averaged over replicates. The last column
#' contains their standard deviation
#' @param maxClust Maximal number of cluster. The minimum is 3
#' @param scaling Either `standardize` (default), `center` or `none`. Standardized 
#' features get mean 0 and standard deviation 1. Centered samples get mean 0. 
#' @param cores The number of threads to be used for parallelisation
#' @return list with the items `ClustInd`: list of clustering objects for each
#' number of clusters, `p` plot object with plots for validity indices,
#' `numclust` optimal cluster number according to "minimum centroid distance"
#' @examples
#' data <- matrix(rnorm(1000), nrow=100)
#' estim_out <- estimClustNum(data, maxClust=10)
#' best_number <- which.max(estim_out[1]) + 2
#' @import limma
#' @import parallel
#' @import stats
#' @import graphics
#' @import grDevices
#' @importFrom shiny getDefaultReactiveDomain incProgress
#' @importFrom matrixStats rowMaxs
#' @export
estimClustNum <- function(dat,
                          maxClust = 25,
                          scaling = "standardize",
                          cores = 1) {
    ClustInd <- matrix(NA, nrow = maxClust - 2, ncol = 6)
    if (is.null(rownames(dat)))
        rownames(dat) <- seq_len(nrow(dat))
    tData <- dat[, seq_len(ncol(dat) - 1)]
    colnames(tData) <- NULL
    
    # define parallelization
    cl <- makeCluster(cores)
    clusterExport(
        cl = cl,
        varlist = c("vsclust_algorithm"),
        envir = environment()
    )
    clusterEvalQ(cl = cl, library(vsclust))
    
    sds <- dat[rownames(tData), ncol(dat)]
    
    #Standardize
    # scale standard deviations by the ones in the actual data to cope for the
    # following standardization
    if (!any(scaling == c("standardize", "center", "none")))
        stop("parameter scaling needs to be standardize, center or none!")
    if ((scaling == "standardize"))
        sds <- sds / rowSds(as.matrix(tData), na.rm = TRUE)
    tData <- t(scale(t(tData), center = (scaling != "none"), scale = (scaling == "standardize")))
    
    multiOut <- lapply(seq(3,maxClust,1), function(x) {
        if (!is.null(getDefaultReactiveDomain())) {
            incProgress(1, detail = paste("Running cluster number", x))
        } else {
            message("Running cluster number ", x)
        }
        clustout <- ClustComp(
            tData,
            NClust = x,
            Sds = sds,
            NSs = 16,
            cl = cl
        )
        c(clustout$indices, sum(rowMaxs(clustout$Bestcl$membership) > 0.5),
          sum(rowMaxs(clustout$Bestcl2$membership) > 0.5))
    })
    
    stopCluster(cl)
    
    for (NClust in seq(3,maxClust,1))
        ClustInd[NClust - 2,] <- multiOut[[NClust - 2]]
    rownames(ClustInd) <- paste0("num_clust_", seq(3,maxClust,1))
    colnames(ClustInd) <-
        c(
            "MinCentroidDist_VSClust",
            "XieBeni_VSClust",
            "MinCentroidDist_FCM",
            "XieBeni_FCM",
            "NumVSClust",
            "NumFCM"
        )
    
    # Output
    ClustInd
}

#' Wrapper for running cluster analysis
#'
#' This function runs the clustering and visualizes the results.
#'
#' @param dat matrix or data frame with feature values for different conditions
#' @param NClust Number of cluster for running the clustering
#' @param proteins vector with additional feature information (default is NULL)
#' to be added to the results
#' @param VSClust boolean. TRUE for running the variance-sensitive clustering.
#' Otherwise, the function will call standard fuzzy c-means clustering
#' @param scaling Either `standardize` (default), `center` or `none`. Standardized 
#' features get mean 0 and standard deviation 1. Centered samples get mean 0.
#' @param cores Number of threads for the parallelization
#' @param verbose Show more information during execution
#' @return list with the items `dat`(the original data), `Bestcl` clustering
#' results (same as from vsclust_algorithm), `p` (plot object with mfuzz plots),
#' `outFileClust`(suitable matrix with complete information) , `ClustInd`
#' (information about being member of any cluster, feature needs on membership
#' values > 0.5)
#' @examples
#' data(iris)
#' data <- cbind(iris[,seq_len(4)],1)
#' clust_out <- runClustWrapper(data, NClust=3, cores=1)
#' clust_out$p
#' @import parallel
#' @import graphics
#' @importFrom grDevices recordPlot
#' @importFrom shiny getDefaultReactiveDomain incProgress
#' @importFrom matrixStats rowMaxs
#' @export
runClustWrapper <-
    function(dat,
             NClust,
             proteins = NULL,
             VSClust = TRUE,
             scaling = "standardize",
             cores,
             verbose = FALSE) {
        tData <- dat[, seq_len(ncol(dat) - 1)]
        sds <- dat[, ncol(dat)]
        
        #Standardize
        # scale standard deviations by the ones in the actual data to cope for the
        # following standardization
        if (!any(scaling == c("standardize", "center", "none")))
            stop("parameter scaling needs to be standardize, center or none!")
        if ((scaling == "standardize"))
            sds <- sds / rowSds(as.matrix(tData), na.rm = TRUE)
        tData <- t(scale(t(tData), center = (scaling != "none"), scale = (scaling == "standardize")))
        if (is.null(rownames(tData))) {
            rownames(tData) <- seq_len(nrow(tData))
        }
        cl <- makeCluster(cores)
        clusterExport(
            cl = cl,
            varlist = c("vsclust_algorithm"),
            envir = environment()
        )
        clusterEvalQ(cl = cl, library(vsclust))
        
        clustout <- ClustComp(
            tData,
            NClust = NClust,
            Sds = sds,
            NSs = 16,
            cl = cl,
            verbose = verbose
        )
        stopCluster(cl)
        
        
        if (VSClust) {
            Bestcl <- clustout$Bestcl
        } else {
            Bestcl <- clustout$Bestcl2
        }
        Bestcl <- SwitchOrder(Bestcl, NClust)
        
        # sorting for membership values (globally)
        Bestcl$cluster <-
            Bestcl$cluster[order(rowMaxs(Bestcl$membership, na.rm = TRUE))]
        Bestcl$membership <-
            Bestcl$membership[order(rowMaxs(Bestcl$membership, na.rm = TRUE)),]
        tData <- tData[names(Bestcl$cluster),]
        
        if (!is.null(getDefaultReactiveDomain()))
            incProgress(0.7, detail = paste("Plotting", NClust))
        
        # graphics.off() ## clean up device
        par(lwd = 0.25)
        oldmar <- par("mar")
        par(mar = c(2, 2, 3, 3), mgp = c(2, 1, 0))
        par(mar = par("mar") / max(1, NClust / 20))
        
        mfuzz.plot(
            tData,
            cl = Bestcl,
            mfrow = c(round(sqrt(NClust)), ceiling(sqrt(NClust))),
            minMem = 0.5,
            colo = "fancy"
        )
        p <- recordPlot()
        # par(lwd=1,mar=oldmar)
        
        colnames(Bestcl$membership) <-
            paste("membership of cluster", colnames(Bestcl$membership))
        outFileClust <- tData
        if (!is.null(proteins)) {
            outFileClust <-
                cbind(outFileClust, names =
                          as.character(proteins[rownames(outFileClust)]))
        }
        
        rownames(Bestcl$centers) <-
            paste("Cluster", rownames(Bestcl$centers))
        ClustInd <-
            as.data.frame(table(Bestcl$cluster[rowMaxs(Bestcl$membership) > 0.5]))
        if (ncol(ClustInd) == 2)
            colnames(ClustInd) <- c("Cluster", "Members")
        else
            ClustInd <-
            cbind(seq_len(max(Bestcl$cluster)), rep(0, max(Bestcl$cluster)))
        
        ## Output
        Out <-
            list(
                dat = tData,
                Bestcl = Bestcl,
                p = p,
                outFileClust = outFileClust,
                ClustInd = ClustInd
            )
        return(Out)
    }

#' Functional Enrichment with STRING
#'
#' @description
#' \code{runFuncEnrich} performs a functional enrichment analysis for each cluster
#' of features (genes/proteins) using an \code{enrichSTRING_API} call (which queries
#' the STRING database). It replaces a previously used approach that relied on
#' \emph{RDAVIDWebService}, and is therefore more up to date.
#'
#' @details
#' The function takes a clustering result (e.g., from \code{vsclust_algorithm} or
#' \code{Bestcl} objects) and, for each cluster:
#' \enumerate{
#'   \item Extracts the member features with membership above 0.5.
#'   \item Optionally replaces their IDs with entries from \code{protnames}.
#'   \item Calls \code{compareCluster} (from \emph{clusterProfiler}) using a
#'         custom \code{enrichSTRING_API} function for the actual STRING-based enrichment.
#' }
#' The resulting \code{compareClusterResult} includes adjusted p-values (FDR),
#' and the top 20 enriched terms are retained if the total set is larger than 20.
#'
#' @param cl A clustering result. Typically either the direct output of
#'        \code{vsclust_algorithm} or a \code{Bestcl} object from \code{ClustComp}
#'        or \code{runClustWrapper}.
#' @param protnames Optional named vector mapping feature IDs (as in \code{cl$cluster})
#'        to more interpretable gene/protein identifiers. If \code{NULL} (the default),
#'        the feature names in \code{cl} themselves are used.
#' @param infosource Currently unused; previously indicated enrichment categories
#'        when using DAVID. Kept for compatibility but is ignored in this version.
#'
#' @return A list with three components:
#' \describe{
#'   \item{\code{fullFuncs}}{The full \code{compareCluster} result from \code{enrichSTRING_API}.}
#'   \item{\code{redFuncs}}{A \code{compareClusterResult} object containing only
#'         the top 20 enriched terms (if more than 20 were detected).}
#'   \item{\code{BHI}}{A numeric value from \code{calcBHI} measuring cluster
#'         heterogeneity.}
#' }
#'
#' @seealso
#' \code{\link{compareCluster}}, \code{\link{enrichSTRING_API}}
#'
#' @examples
#' \dontrun{
#'   # Suppose 'mycl' is a clustering result from vsclust_algorithm,
#'   # and we have a named vector 'myProtnames' that maps feature IDs to gene symbols:
#'
#'   res <- runFuncEnrich(mycl, protnames = myProtnames)
#'   res$fullFuncs  # The full compareCluster output
#'   res$redFuncs   # The top 20 enriched terms
#'   res$BHI        # The numeric BHI value
#' }
#'
#' @importFrom clusterProfiler compareCluster
#'
#' @export
runFuncEnrich <-
    function(cl, protnames = NULL, infosource) {
        Accs <- list()
        for (c in seq_len(max(cl$cluster))) {
            cname <- paste("Cluster", c, sep = "_")
            Accs[[cname]] <-
                names(which(cl$cluster == c & rowMaxs(cl$membership) > 0.5))
            
            Accs[[cname]] <- Accs[[cname]][Accs[[cname]] != ""]
            if (length(Accs[[cname]]) > 0) {
                if (!is.null(protnames)) {
                    Accs[[cname]] <- as.character(protnames[Accs[[cname]]])
                }
                
                Accs[[cname]] <- sub("-[0-9]", "", Accs[[cname]])
            }
        }
        # TODO? add extraction of multiple accession numbers
        Accs <- lapply(Accs, function(x)
            unique(ifelse(is.na(x), "B3", x)))
        Accs <- Accs[lapply(Accs, length) > 0]
        x <- NULL
        try(x <-
                clusterProfiler::compareCluster(
                    Accs,
                    fun = "enrichSTRING_API",
                    category = unlist(infosource)
                ))
        validate(need(!is.null(x), "No result. IDs might not be matching or you do not have any significant enrichments"))
        if (!is.null(getDefaultReactiveDomain()))
            incProgress(0.7, detail = "received")
        message("got data from STRINGdb\n")
        x@compareClusterResult <- cbind(x@compareClusterResult,
                                        log10padval =
                                            log10(x@compareClusterResult$p.adjust))
        y <-
            new("compareClusterResult",
                compareClusterResult = x@compareClusterResult)
        if (length(unique(y@compareClusterResult$ID)) > 20) {
            message("Reducing number of STRINGdb results\n")
            y@compareClusterResult <- y@compareClusterResult[order(y@compareClusterResult$p.adjust)[seq_len(20)],]
            
            y@compareClusterResult$Cluster <-
                as.character(y@compareClusterResult$Cluster)
        }
        
        BHI <- calcBHI(Accs, x)
        return(list(
            fullFuncs = x,
            redFuncs = y,
            BHI = BHI
        ))
        
    }




#' Enrichment Analysis via STRING REST API
#'
#' @description
#' \code{enrichSTRING_API} performs a functional enrichment analysis by sending a
#' gene list to STRING's TSV API endpoint, retrieving enrichment results for one or
#' more categories (e.g., "KEGG"), and building a \code{clusterProfiler}-style result
#' object. It is intended as a lightweight replacement for older web-service–based
#' methods like \emph{RDAVIDWebService}.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Accepts a vector of \code{genes} recognized by STRING (e.g., Ensembl,
#'         UniProt, or commonly used gene symbols for a given species).
#'   \item Sends a \emph{POST} request to \url{https://string-db.org/api/tsv/enrichment} with
#'         those identifiers (and optionally \code{species}).
#'   \item Parses the returned \emph{TSV}-formatted enrichment data, which typically
#'         includes columns for \code{term}, \code{description}, \code{p_value}, and more.
#'   \item Optionally filters by \code{category} (e.g., "KEGG", "Process", etc.), applies
#'         \code{BH} multiple-testing correction, and removes terms above the
#'         \code{adjpvalueCutoff}.
#'   \item Returns an object of class \code{enrichResult} if \emph{clusterProfiler} is
#'         installed; otherwise, it simply returns the filtered data frame.
#' }
#'
#' Note that this function does \emph{not} allow you to explicitly supply a custom
#' background gene set. The STRING API by default uses the entire known set of genes
#' or proteins for the specified species as the background. Also, the total number of
#' background genes is not reported by STRING; only how many of them map to each term.
#'
#' @param genes Character vector of gene or protein identifiers recognized by STRING.
#' @param species A single numeric or string specifying the NCBI taxonomy ID
#'        (e.g., \code{9606} for human). If \code{"none"} (default), no species is set in
#'        the request, and STRING will attempt to auto-detect or may fail.
#' @param category One or more enrichment categories from STRING (e.g., \code{"KEGG"},
#'        \code{"Process"}, etc.). Defaults to \code{"KEGG"} if unspecified.
#' @param adjpvalueCutoff Numeric cutoff for the BH-adjusted p-value (default 0.05).
#' @param verbose Logical indicating whether to print diagnostic messages (default \code{FALSE}).
#'
#' @return If any terms pass the cutoff, an \code{enrichResult} object
#'   (from \emph{clusterProfiler}) is returned. Otherwise, if \emph{clusterProfiler}
#'   is not installed or no terms pass filtering, the function either returns a
#'   data.frame (if installed but no terms pass) or \code{NULL}.
#'
#' @seealso
#' \href{https://string-db.org/help/api/}{STRING API documentation},
#' \code{\link{clusterProfiler}}.
#'
#' @examples
#'   library(httr)
#'
#'   # A small gene set:
#'   gene_set <- c("TP53","BRCA1","BRCA2","EGFR")
#'
#'   # Perform enrichment on KEGG terms for human (9606):
#'   enr <- enrichSTRING_API(
#'             genes       = gene_set,
#'             species     = 9606,
#'             category    = "KEGG",
#'             adjpvalueCutoff = 0.05,
#'             verbose     = TRUE
#'          )
#'
#'   if (!is.null(enr)) {
#'     # If clusterProfiler is installed and some terms pass filtering, check results:
#'     head(enr@result)
#'   }
#'   
#' @import httr 
#' @import DOSE
#'
#' @export
enrichSTRING_API <- function(genes,
                             species        = "none",
                             category       = "KEGG",
                             adjpvalueCutoff   = 0.05,
                             verbose        = FALSE
                             ) {
    # 1) Check for required packages
    if (!requireNamespace("httr", quietly=TRUE)) {
        stop("Package 'httr' is required but not installed.")
    }
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        # needed if you want to build an 'enrichResult' object
        warning("Package 'clusterProfiler' not installed; building 'enrichResult' might fail.")
    }

    # 2) Build POST request body
    # According to the API docs, identifiers are line-separated. We'll use '\r' or '\n'
    identifiers_str <- paste(genes, collapse="%0d")

    # Build a list for the POST fields
    body_list <- list(
        identifiers     = identifiers_str,
        caller_identity = "VSClust_enrichSTRING_API"  # customize as you like
    )
    print(body_list)
    if (species != "none") {
        body_list$species <- species
    }    
        
    # The STRING endpoint:
    url <- "https://string-db.org/api/tsv/enrichment"

    # 3) Make the POST request
    if (verbose) message("Contacting STRING for functional enrichment via POST...")
    resp <- httr::POST(url = url, body = body_list, encode = "form")
    
    # 4) Check status code
    code <- httr::status_code(resp)
    if (code != 200) {
        # Try to retrieve error message text
        err_msg <- httr::content(resp, as = "text", encoding = "UTF-8")
        stop("STRING API request failed. Status code: ", code, "\n", err_msg)
    }

    # 5) Parse the TSV content from the response
    resp_text <- httr::content(resp, as="text", encoding="UTF-8")
    # If there's an error from STRING, it might say something like "ERROR ..."
    if (grepl("ERROR", resp_text, ignore.case=TRUE)) {
        warning("STRING returned an error or no enrichment found:\n", resp_text)
        return(NULL)
    }
    
    # Read the text as TSV
    # The columns are typically: 
    # category,term,description,input_number,background_number,...
    # p_value,fdr,preferred_names, ...
    # We'll attempt to parse them with read.table:
    df <- tryCatch({
        utils::read.table(textConnection(resp_text),
                          sep="\t", header=TRUE, quote="", comment.char="",
                          stringsAsFactors=FALSE)
    }, error=function(e) {
        stop("Failed to parse STRING API response as TSV:\n", e$message)
    })
    
    # 6) Check if data frame is empty or if required columns are missing
    needed_cols <- c("term","description","p_value","fdr","category",
                     "number_of_genes","number_of_genes_in_background")
    if (!all(needed_cols %in% colnames(df))) {
        warning("Some expected columns are missing from the STRING response. Possibly no results.")
        return(NULL)
    }
    if (nrow(df) == 0) {
        warning("STRING returned an empty table (no enrichment).")
        return(NULL)
    }
    
    # 7) Filter for requested category (if you provided multiple, we can subset)
    if (!is.null(category) && length(category) >= 1) {
        df <- df[df$category %in% category, , drop=FALSE]
    }
    if (nrow(df) == 0) {
        warning("No terms remain after filtering by category.")
        return(NULL)
    }

    # 8) Build the result table, similar to an "enrichResult"
    Over <- data.frame(
        ID          = df$term,           # e.g. "GO:0006260"
        Description = df$description,    # e.g. "DNA replication"
        pvalue      = as.numeric(df$p_value),
        fdr         = as.numeric(df$fdr),
        Count   = as.numeric(df$number_of_genes),  # how many from our input
        bgCount     = as.numeric(df$number_of_genes_in_background),# how many in STRING's background
        Category    = df$category,       # e.g. "Process", "KEGG", ...
        geneID      = df$inputGenes,
        stringsAsFactors = FALSE,
        GeneRatio = 1  # we don't have the background count, so assume 1:1
    )
    
    # 9) We'll create p.adjust from pvalue ourselves, ignoring the 'fdr' from STRING if we want
    Over$p.adjust <- stats::p.adjust(Over$pvalue, method="BH")
    
    # 11) Filter by SSize, maxGSSize, pvalueCutoff, etc.
    keep_padjust   <- Over$p.adjust <= adjpvalueCutoff
    keep_idx <- keep_padjust
    
    Over <- Over[keep_idx, , drop=FALSE]
    
    if (nrow(Over) == 0) {
        warning("No terms passed the final filtering steps.")
        return(NULL)
    }
    
    # 11) Return the final data.frame
    #     Or you could build an S4 'enrichResult' object if you want clusterProfiler-like usage.
    new("enrichResult", result = Over, pvalueCutoff = adjpvalueCutoff, 
        pAdjustMethod = "BH", organism = as.character(unique(df$ncbiTaxonId)), ontology = "STRINGout", 
        keytype = "varied")
}


#' Run VSClust as Shiny app
#'
#' You will get the full functionality of the VSClust workflow with multiple
#' visualizations and downloads
#'
#' @return The shiny app should open in a browser or in RStudio.
#' @examples
#' \dontrun{
#' runVSClustApp()}
#' @export
#' @references
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering
#' of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi:
#' 10.1093/bioinformatics/bty224. PMID: 29635359.
#'
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and
#' the Quantitative Analysis of Protein Complexes. Methods Mol Biol.
#' 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#'
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters
#' for fuzzy c-means cluster analysis. Bioinformatics. 2010
#' Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29.
#' PMID: 20880957.
runVSClustApp <- function() {
    shiny::runApp(system.file("shiny/", package = "vsclust"))
}
