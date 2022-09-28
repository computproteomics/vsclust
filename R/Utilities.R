validate <- shiny::validate

# extend to 702 cases:
LETTERS702 <-
  c(LETTERS, vapply(LETTERS, function(x)
    paste0(x, LETTERS), character(26)))


# Some mathematical functions
#' @import stats
erf <- function(x)
  2 * pnorm(x * sqrt(2)) - 1
#' @import stats
erf.inv <- function(x)
  qnorm((x + 1) / 2) / sqrt(2)


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

#' arrange cluster member numbers from largest to smallest
#'
#' @param Bestcl fclust object
#' @param NClust Number of clusters
#' @importFrom matrixStats rowMaxs
#' @return fclust object with reorder clusters
#' @example 
#' @examples
#' # Generate some random data
#' data <- matrix(rnorm(seq_len(1000)), nrow=100)
#' # Run clustering
#' clres <- vsclust_algorithm(data, centers=5, m=1.5)
#' clres <- SwitchOrder(clres, 5)
#' @export
SwitchOrder <- function(Bestcl, NClust) {
  switching <-
    as.numeric(names(sort(table(Bestcl$cluster[rowMaxs(Bestcl$membership) >
                                                 0.5]),
                          decreasing = TRUE)))
  if (length(switching) < NClust)
    switching <-
      c(switching, which(!((seq_len(
        NClust
      )) %in% switching)))
  switching2 <- seq_len(NClust)
  names(switching2) <- switching
  tBest <- Bestcl
  tBest$centers <- Bestcl$centers[switching,]
  rownames(tBest$centers) <- seq_len(NClust)
  tBest$size <- Bestcl$size[switching]
  tBest$cluster <- switching2[as.character(Bestcl$cluster)]
  names(tBest$cluster) <- names(Bestcl$cluster)
  tBest$membership <- Bestcl$membership[, switching]
  colnames(tBest$membership) <- seq_len(NClust)
  tBest
}

#' Calculate "biological homogeneity index"
#'
#' This index is providing a number for the enriched GO terms and pathways to
#' assess the biological content within a set of genes or proteins.
#' The calculation is according to Datta, S. & Datta, S. Methods for evaluating
#' clustering algorithms for gene expression data using a reference set of
#' functional classes. BMC bioinformatics 7, 397 (2006).
#' @param Accs list containing gene or protein IDs, such as UniProt accession
#' names
#' @param gos object from ClusterProfiler
#' @return Biological Homogeneity Index
#' @examples
#' # Run enrichment analysis
#' library(clusterProfiler)
#' data(gcSample)
#' xx <- compareCluster(gcSample, fun="enrichKEGG",
#'                      organism="hsa", pvalueCutoff=0.05)
#' # Generate random list from gcSample
#' rand_ids <- lapply(gcSample, function(x) sample(unlist(gcSample), 200))
#' calcBHI(rand_ids, xx)
#'
#' @import stats
#' @export
#' @references
#' Datta, S. & Datta, S. Methods for evaluating
#' clustering algorithms for gene expression data using a reference set of
#' functional classes. BMC bioinformatics 7, 397 (2006).
#'
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering
#' of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972.
#' doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#'
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and
#' the Quantitative Analysis of Protein Complexes. Methods Mol Biol.
#' 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#'
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters
#' for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8.
#' doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.

calcBHI <- function(Accs, gos) {
  ## enrichment might not yield all GO terms! This could lead to problems
  BHI <- sumcomb <- vector("integer", length(Accs))
  names(BHI) <- names(Accs)
  goData <- as.data.frame(gos@compareClusterResult)
  for (i in names(gos@geneClusters)) {
    genes <- Accs[[i]]
    ispair <-
      matrix(FALSE,
             length(genes),
             length(genes),
             dimnames = list(rows = genes, cols = genes))
    sumcomb[i] <- choose(length(genes), 2)
    clgroup <- goData[goData$Cluster == i, "geneID"]
    for (j in seq_len(length(clgroup))) {
      tgenes <-
        genes[na.omit(match(
          unlist(strsplit(clgroup[j], "/"), use.names = FALSE),
          as.character(genes)
        ))]
      ltgenes <- length(tgenes)
      if (ltgenes > 1) {
        for (i1 in tgenes[seq_len(ltgenes - 1)]) {
          ttgene <- tgenes[(which(i1 == tgenes) + 1):ltgenes]
          ispair[i1, ttgene] <- ispair[ttgene, i1] <- TRUE
        }
        
      }
      # print(tgenes)
    }
    BHI[i] <- BHI[i] + sum(ispair) / 2
  }
  # Rprof()
  # print(summaryRprof(tmp))
  sum(BHI) / sum(sumcomb)
}

#' Calculate mean over replicates
#'
#' Simple method to calculate the means for each feature across its replicates
#'
#' @param data Matrix of data frame with numerical values. Columns corresponds
#' to samples
#' @param NumReps Number of replicates per experimental condition
#' @param NumCond Number of different experimental conditions
#' @return Matrix of data frame with averaged values over replicates for each
#' conditions
#' @examples
#' data <- matrix(rnorm(1000), nrow=100)
#' av_data <- averageCond(data, NumCond=2, NumReps=5)
#' @export
averageCond <- function(data, NumReps, NumCond) {
  # Calculates means over replicates
  tdat <-
    rowMeans(data[, seq(1, NumReps * NumCond, NumCond)], na.rm = TRUE)
  for (i in 2:NumCond) {
    tdat <-
      cbind(tdat, rowMeans(data[, seq(i, NumReps * NumCond, NumCond)], na.rm =
                             TRUE))
  }
  colnames(tdat) <-
    paste("Mean of log ", LETTERS702[seq_len(NumCond)], sep = "")
  tdat
  
}

# Function to add columns with NAs for balancing number of replicates
# Also rearranges the columns to replicate groups
# Create data set
balanceData <- function(dat, coldat) {
  tdat <- NULL
  NumReps <- max(table(coldat))
  NumCond <- length(unique(coldat))
  for (i in unique(coldat)) {
    tdat <- cbind(tdat, dat[, names(coldat == i)])
    # add empty columns to achieve the same number of replicates per condition
    if (length(coldat == i) < NumReps) {
      tdat <- cbind(tdat, matrix(
        NA,
        nrow = nrow(tdat),
        ncol = NumReps - length(coldat == i)
      ))
    }
                                                           }
  # Rearrange columns
  tdat[, rep(seq_len(NumCond)-1, NumReps) * NumReps + rep(seq_len(NumReps),
                                                       each = NumCond)]
}
