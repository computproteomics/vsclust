#' 
#' Determine individual fuzzifier values
#'
#' This function calculated the values of the fuzzifier from a) the dimensions
#' of the considered data set and b)
#' from the individual feature standard deviations.
#'
#' @param dims vector of two integers containing the dimensions of the data 
#' matrix for the clustering
#' @param NClust Number of cluster for running vsclust on (does no influence the 
#' calculation of `mm`)
#' @param Sds individual standard deviations, set to 1 if not available
#' @return list of `m`: individual fuzzifiers, `mm`: standard fuzzifier for fcm 
#' clustering when not using vsclust algorithm
#' @examples
#' # Generate some random data
#' data <- matrix(rnorm(seq_len(1000)), nrow=100)
#' # Estimate fuzzifiers
#' fuzz_out <- determine_fuzz(dim(data), 1)
#' # Run clustering
#' clres <- vsclust_algorithm(data, centers=5, m=fuzz_out$mm)
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
#' for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;
#' 26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. 
#' PMID: 20880957.
determine_fuzz <- function(dims, NClust, Sds = 1) {
  D <- dims[2]
  d <- sqrt(D / 2)
  mm <-
    1 + (1418 / dims[1] + 22.05) * dims[2] ^ (-2) + (12.33 / dims[1] + 0.243) *
    dims[2] ^ (-0.0406 * log(dims[1]) - 0.1134)
  
  ### d_i and d_t
  difunc <-
    function(c, D) {
      x <- 0:c
      sum(choose(c, x) / (x * D + 1) * (-1) ^ x)
    }
  
  di <- difunc(NClust, D)  / sqrt(pi) * exp(lgamma(D / 2 + 1) * (1 / D))
  dt <- (NClust) ^ (-1 / D)
  p <- dnorm(di, 0, Sds) * (1 - dnorm(dt, 0, Sds)) ^ (NClust - 1)
  
  m <- mm + p * mm * (D / 3 - 1)
  m[m == Inf] <- 0
  m[m == 0] <- NA
  m[is.na(m)] <- mm * 10
  
  ## If m for highest Sd is mm then all = mm
  if (m[which.max(Sds)] == mm)
    m[seq_len(length(m))] <- mm
  
  return(list (m = m, mm = mm))
  
}

#' Run the vsclust clustering algorithm
#'
#' This function calls the c++ implementation of the vsclust algorithm, being an 
#' extension of fuzzy c-means clustering with additional variance control and 
#' capability to run on data with missing values
#'
#' @param x a numeric data matrix
#' @param centers Either numeric for number of clusters or numeric matrix with 
#' center coordinates
#' @param iter.max Numeric for maximum number of iterations
#' @param verbose Verbose information
#' @param dist Distance to use for the calculation. We prefer "euclidean" 
#' (default)
#' @param m Fuzzifier value: numeric or vector of length equal to number of rows 
#' of x
#' @param rate.par (experimental) numeric value for punishing missing values
#' @param weights numeric or vector of length equal to number of rows of x
#' @param control list with arguments to vsclust algorithms (now only cutoff for 
#' relative tolerance: reltol)
#' @return list with details about clustering having the objects `centers` 
#' (positions of centroids), `size` (feature number per cluster), `cluster` 
#' (nearest cluster of each feature), `membership` matrix of membership values, 
#' `iter` (number of carried out iterations), `withinerror` (final error from 
#' optimization), `call`(call of function)
#' @examples
#' #' # Generate some random data
#' data <- matrix(rnorm(seq_len(1000)), nrow=100)
#' # Run clustering
#' clres <- vsclust_algorithm(data, centers=5, m=1.5)
#' head(clres$membership)
#' @export
#' @useDynLib
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
vsclust_algorithm <-
  function(x,
           centers,
           iter.max = 100,
           verbose = FALSE,
           dist = "euclidean",
           m = 2,
           rate.par = NULL,
           weights = 1,
           control = list())
  {
    x <- as.matrix(x)
    xrows <- nrow(x)
    xcols <- ncol(x)
    
    if (missing(centers))
      stop("Argument 'centers' must be a number or a matrix.")
    
    dist <- pmatch(dist, c("euclidean", "manhattan"))
    if (is.na(dist))
      stop("invalid distance")
    if (dist == -1)
      stop("ambiguous distance")
    
    if (length(centers) == 1) {
      ncenters <- centers
      centers <- x[sample(seq_len(xrows), ncenters), , drop = FALSE]
      centers[is.na(centers)] <- 0
      if (any(duplicated(centers))) {
        cn <- unique(x)
        mm <- nrow(cn)
        if (mm < ncenters)
          stop("More cluster centers than distinct data points.")
        centers <- cn[sample(seq_len(mm), ncenters), , drop = FALSE]
      }
    } else {
      centers <- as.matrix(centers)
      if (any(duplicated(centers)))
        stop("Initial centers are not distinct.")
      cn <- NULL
      ncenters <- nrow(centers)
      if (xrows < ncenters)
        stop("More cluster centers than data points.")
    }
    
    if (xcols != ncol(centers))
      stop("Must have same number of columns in 'x' and 'centers'.")
    
    if (iter.max < 1)
      stop("Argument 'iter.max' must be positive.")
    
    if (missing(rate.par)) {
      rate.par <- 0
    }
    
    reltol <- control$reltol
    if (is.null(reltol))
      reltol <- sqrt(.Machine$double.eps)
    if (reltol <= 0)
      stop("Control parameter 'reltol' must be positive.")
    
    if (any(weights < 0))
      stop("Argument 'weights' has negative elements.")
    if (!any(weights > 0))
      stop("Argument 'weights' has no positive elements.")
    weights <- rep(weights, length = xrows)
    weights <- weights / sum(weights)
    
    # if length of fuzzifiers is lower than number of features, repeat the 
    # pattern until end.
    # Also counts for single fuzzifier
    m <- rep(m, length = xrows)
    
    initcenters <- centers
    pos <- as.factor(seq_len(ncenters))
    rownames(centers) <- pos
    
    u <- matrix(0.0,  nrow = xrows, ncol = ncenters)
    iter <- c(0L)
    val <-
      c_plusplus_means(x,
                       centers,
                       weights,
                       m,
                       dist - 1,
                       iter.max,
                       reltol,
                       verbose,
                       u ,
                       1,
                       iter,
                       NA,
                       rate.par)
    # put modified values in retval
    retval <-
      list(
        x = x,
        xrows = xrows,
        xcols = xcols,
        centers = centers,
        ncenters = ncenters,
        m = m,
        dist = dist - 1,
        iter.max = iter.max,
        reltol = reltol,
        verbose = verbose,
        rate.par = rate.par,
        u = u,
        ermin = val,
        iter = iter
      )
    
    centers <- matrix(retval$centers,
                      ncol = xcols,
                      dimnames = list(seq_len(ncenters),
                                      colnames(initcenters)))
    u <- matrix(retval$u,
                ncol = ncenters,
                dimnames = list(rownames(x), seq_len(ncenters)))
    # u <- u[order(perm), ]
    iter <- retval$iter - 1
    withinerror <- retval$ermin
    
    cluster <- apply(u, 1, which.max)
    clustersize <- as.integer(table(cluster))
    
    retval <- list(
      centers = centers,
      size = clustersize,
      cluster = cluster,
      membership = u,
      iter = iter,
      withinerror = withinerror,
      call = match.call()
    )
    class(retval) <- c("fclust")
    return(retval)
  }

#' Function to run clustering with automatic fuzzifier settings (might become 
#' obsolete)
#'
#' Run original fuzzy c-means and vsclust for a number of clusters and the given 
#' data set including data pre-processing and
#' automatic setting of the data-dependent parameters like the lower limit of 
#' the fuzzifier.
#'
#' @param dat a numeric data matrix
#' @param NSs number of clusterings runs with different random seeds
#' @param NClust Number of clusters
#' @param Sds Standard deviation of features (either vector of the same length 
#' as features numbers in matrix or single value)
#' @param cl object of class `cluster` or `SOCKcluster` to specify environment 
#' for parallelization
#' @return List containing the objects
#' @return `indices` containing minimum centroid distance and Xie-Beni index for 
#' both clustering methods
#' @return `Bestcl` optimal vsclust results (variance-sensitive fcm clustering)
#' @return `Bestcl2` optimal fuzzy c-means restults
#' @return `m` vector of individual fuzzifer values per feature
#' @return `withinerror` final optimization score for vsclust
#' @return `withinerror2` final optimization score for fuzzy c-means clustering
#' @examples
#' #' # Generate some random data
#' data <- matrix(rnorm(seq_len(1000)), nrow=100)
#' # Run clustering
#' cl <- parallel::makePSOCKcluster(1, nnodes=1)
#' ClustCompOut <- ClustComp(data, cl=cl, NClust=6, Sds=1)
#' barplot(ClustCompOut$indices)
#' @import parallel
#' @import stats
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
ClustComp <-
  function(dat,
           NSs = 10,
           NClust = NClust,
           Sds = Sds,
           cl = parallel::makePSOCKcluster(1)) {
    fuzz_out <- determine_fuzz(dim(dat), NClust, Sds)
    m <- fuzz_out$m
    mm <- fuzz_out$mm
    
    clusterExport(
      cl = cl,
      varlist = c("dat", "NClust", "m", "mm"),
      envir = environment()
    )
    cls <-
      parLapply(cl, seq_len(NSs), function(x)
        vsclust_algorithm(
          dat,
          NClust,
          m = m,
          verbose = FALSE,
          iter.max =
            1000
        ))
    # cls <- lapply(seq_len(NSs), function(x) vsclust_algorithm(tData,NClust,
    # m=m, verbose=FALSE,iter.max=1000))  #print(cls[[1]])
    Bestcl <- cls[[which.min(lapply(cls, function(x)
      x$withinerror))]]
    cls <-
      parLapply(cl, seq_len(NSs), function(x)
        vsclust_algorithm(
          dat,
          NClust,
          m = mm,
          verbose = FALSE,
          iter.max =
            1000
        ))
    Bestcl2 <- cls[[which.min(lapply(cls, function(x)
      x$withinerror))]]
    # stopCluster(cl)
    
    # return validation indices
    list(
      indices = c(
        min(dist(Bestcl$centers)),
        cvalidate.xiebeni(Bestcl, mm),
        min(dist(Bestcl2$centers)),
        cvalidate.xiebeni(Bestcl2, mm)
      ),
      Bestcl = Bestcl,
      Bestcl2 = Bestcl2,
      m = m,
      withinerror = Bestcl$withinerror,
      withinerror2 = Bestcl2$withinerror
    )
  }
