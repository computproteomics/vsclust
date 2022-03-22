vsclust_algorithm <-
  function(x, centers, iter.max = 100, verbose = FALSE,
           dist = "euclidean", m = 2,
           rate.par = NULL, weights = 1, control = list())
  {
    x <- as.matrix(x)
    xrows <- nrow(x)
    xcols <- ncol(x)
    
    if(missing(centers))
      stop("Argument 'centers' must be a number or a matrix.")
    
    dist <- pmatch(dist, c("euclidean", "manhattan"))
    if(is.na(dist)) 
      stop("invalid distance")
    if(dist == -1) 
      stop("ambiguous distance")
    
    if(length(centers) == 1) {
      ncenters <- centers
      centers <- x[sample(1 : xrows, ncenters), , drop = FALSE]
      centers[is.na(centers)] <- 0
      if(any(duplicated(centers))) {
        cn <- unique(x)
        mm <- nrow(cn)
        if(mm < ncenters) 
          stop("More cluster centers than distinct data points.")
        centers <- cn[sample(1 : mm, ncenters), , drop = FALSE]
      }
    } else {
      centers <- as.matrix(centers)
      if(any(duplicated(centers))) 
        stop("Initial centers are not distinct.")
      cn <- NULL
      ncenters <- nrow(centers)
      if (xrows < ncenters)
        stop("More cluster centers than data points.")
    }
    
    if(xcols != ncol(centers))
      stop("Must have same number of columns in 'x' and 'centers'.")
    
    if(iter.max < 1) 
      stop("Argument 'iter.max' must be positive.")
    
    if(missing(rate.par)) {
      rate.par <- 0
    }
    
    reltol <- control$reltol
    if(is.null(reltol))
      reltol <- sqrt(.Machine$double.eps)
    if(reltol <= 0)
      stop("Control parameter 'reltol' must be positive.")
    
    if(any(weights < 0))
      stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
      stop("Argument 'weights' has no positive elements.")
    weights <- rep(weights, length = xrows)
    weights <- weights / sum(weights)
    
    # if length of fuzzifiers is lower than number of features, repeat the pattern until end. 
    # Also counts for single fuzzifier
    m <- rep(m, length = xrows)
    
    initcenters <- centers
    pos <- as.factor(1 : ncenters)
    rownames(centers) <- pos
    
    u <- matrix(0.0,  nrow = xrows, ncol = ncenters)
    iter <- c(0L)
    val <- c_plusplus_means(x, centers, weights, m, dist-1, iter.max, reltol, verbose, u , 1, iter, NA, rate.par)
    # put modified values in retval
    retval <- list(x=x, xrows = xrows, xcols = xcols, centers = centers,ncenters=ncenters, m = m, dist = dist -1,
                   iter.max = iter.max, reltol = reltol, verbose = verbose, rate.par = rate.par, u = u, ermin = val, iter = iter)
    
    centers <- matrix(retval$centers, ncol = xcols,
                      dimnames = list(1 : ncenters,
                                      colnames(initcenters)))
    u <- matrix(retval$u, ncol = ncenters,
                dimnames = list(rownames(x), 1 : ncenters))
    # u <- u[order(perm), ]
    iter <- retval$iter - 1
    withinerror <- retval$ermin
    
    cluster <- apply(u, 1, which.max)
    clustersize <- as.integer(table(cluster))
    
    retval <- list(centers = centers, size = clustersize,
                   cluster = cluster, membership = u, iter = iter,
                   withinerror = withinerror, call = match.call())
    class(retval) <- c("fclust")
    return(retval)
  }

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
