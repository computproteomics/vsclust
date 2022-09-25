#' @section VSClust functions:
#' Functions for running VSClust analysis
#'
#' @docType package
#' @name vsclust
#' @useDynLib vsclust
NULL
#> NULL
#' Wrapper for statistical analysis
#'
#' Prepare data for running vsclust clustering.
#' This includes visualization running the functions for the principal component
#' analysis and its visualization, statistical testing with LIMMA, as well as
#' scaling and filtering of missing values
#' @param dat matrix or data frame of numerical data. Columns are samples.
#' Replicates are grouped (i.e. A1, B1, C1, A2, B2, C2) when letters denote
#' conditions and numbers the replicates. In case of `isStat=F`, you need a
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
#' @param SummarizedExperiment object
#' @param Sample in SummarizedExperiment object
#' @param Column in colData for extracting replicates
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
#' data(miniACC)
#' 
#' stats <- PrepareSEForVSClust(miniACC, coldatname="COC", isStat=TRUE)
#'
#' @import stats
#' @import MultiAssayExperiment
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
#' @param cores The number of threads to be used for parallelisation
#' @return list with the items `ClustInd`: list of clustering objects for each
#' number of clusters, `p` plot object with plots for validity indices,
#' `numclust` optimal cluster number according to "minimum centroid distance"
#' @examples
#' data <- matrix(rnorm(1000), nrow=100)
#' estim_out <- estimClustNum(data, maxClust=10)
#' best_number <- max(estim_out[1])
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
  
  # Standardise
  sds <- dat[rownames(tData), ncol(dat)]
  # scale standard deviations by the ones in the actual data to cope for the
  # following standardization
  sds <- sds / (rowSds(as.matrix(tData), na.rm = TRUE))
  tData <- t(scale(t(tData)))
  
  multiOut <- lapply(3:maxClust, function(x) {
    if (!is.null(getDefaultReactiveDomain())) {
      incProgress(1, detail = paste("Running cluster number", x))
    } else {
      message("Running cluster number", x)
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
  
  for (NClust in 3:maxClust)
    ClustInd[NClust - 2,] <- multiOut[[NClust - 2]]
  rownames(ClustInd) <- paste0("num_clust_", 3:maxClust)
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
#' @param cores Number of threads for the parallelization
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
           cores) {
    tData <- dat[, seq_len(ncol(dat) - 1)]
    sds <- dat[, ncol(dat)]
    
    #Standardize
    # scale standard deviations by the ones in the actual data to cope for the
    # following standardization
    sds <- sds / rowSds(as.matrix(tData), na.rm = TRUE)
    tData <- t(scale(t(tData)))
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
      cl = cl
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
      min.mem = 0.5,
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

#### The manual of the following function was removed to avoid calling
#### RDAVIDWEBSERVICE
# Wrapper for functional enrichment
#
# The functional analysis uses the libarary RDAVIDWebService and thus might
# become obsolete as that library is not supported anymore
# The user can select different ID types and different enrichment categories
# like GO terms and pathways.
# Allowed ID types:
# "AFFYMETRIX_3PRIME_IVT_ID",
# "AFFYMETRIX_EXON_GENE_ID", "AGILENT_CHIP_ID",
# "AGILENT_ID", "AGILENT_OLIGO_ID", "APHIDBASE_ID", "BEEBASE_ID",
# "BEETLEBASE_ID", "BGD_ID", "CGNC_ID", "CRYPTODB_ID", "DICTYBASE_ID",
# "ENSEMBL_GENE_ID",
# "ENSEMBL_TRANSCRIPT_ID", "ENTREZ_GENE_ID", "GENOMIC_GI_ACCESSION",
# "FLYBASE_GENE_ID", "GENBANK_ACCESSION",
# "GENPEPT_ACCESSION", "LOCUS_TAG", "ILLUMINA_ID", "MGI_ID", "MIRBASE_ID",
# "OFFICIAL_GENE_SYMBOL", "PFAM_ID", "PIR_ID", "PROTEIN_GI_ACCESSION",
# "MRNA_GI_ACCESSION",
# "REFSEQ_GENOMIC", "REFSEQ_MRNA", "REFSEQ_PROTEIN", "REFSEQ_RNA",
# "RGD_ID", "SGD_ID", "TAIR_ID", "UCSC_GENE_ID", "UNIGENE",
# "UNIPROT_ACCESSION", "UNIPROT_ID", "UNIREF100_ID", "WORMBASE_GENE_ID",
# "WORMPEP_ID", "ZFIN_ID"
# Allowed enrichment categories:
# "GOTERM_MF_ALL","GOTERM_BP_ALL",
# "GOTERM_CC_ALL","GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"),
# "KEGG"="KEGG_PATHWAY","PANTHER_PATHWAY","REACTOME_PATHWAY","BBID","BIOCARTA",
# "DIP","MINT","INTACT","BIOGRID_INTERACTION","GAD_DISEASE","GAD_DISEASE_CLASS",
# "OMIM_DISEASE","INTERPRO","PROSITE","PFAM","SMART","PRODOM","PIR_SUPERFAMILY"
#
# param cl clustering results (either directly from vsclust_algorithm or as
# `Bestcl` object from ClustComp or runClustWrapper)
# param protnames vector providing the corresponding gene/protein names of the
# features (set to NULL for directly using the feature names (default))
# param idtypes type of IDs for features given by genes/proteins (generic gene
# names are not working)
# param infosource Type of gene annotation (e.g. KEGG_PATHWAY)
# return plot object to be able to pass the figures to e.g. shiny
# @export
runFuncEnrich <-
  function(cl, protnames = NULL, idtypes, infosource) {
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
          compareCluster(
            Accs,
            fun = "enrichDAVID",
            annotation = infosource,
            idType = idtypes,
            david.user = "veits@bmb.sdu.dk"
          ))
    validate(need(!is.null(x), "No result. Wrong ID type?"))
    if (!is.null(getDefaultReactiveDomain()))
      incProgress(0.7, detail = "received")
    message("got data from DAVID\n")
    x@compareClusterResult <- cbind(x@compareClusterResult,
                                    log10padval =
                                      log10(x@compareClusterResult$p.adjust))
    y <-
      new("compareClusterResult",
          compareClusterResult = x@compareClusterResult)
    if (length(unique(y@compareClusterResult$ID)) > 20) {
      message("Reducing number of DAVID results\n")
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

#' Run VSClust as Shiny app
#'
#' You will get the full functionality of the VSClust workflow with multiple
#' visualizations and downloads
#'
#' @return The shiny app should open in a browser or in RStudio.
#' @examples
#' \donttest{
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