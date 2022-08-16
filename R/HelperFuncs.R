#' @section VSClust main functions:
#' The main function for the clustering
#'
#' @docType package
#' @name vsclust
#' @useDynLib vsclust
NULL
#> NULL

########### functions for VSClust


validate <- shiny::validate

# extend to 702 cases:
LETTERS702 <- c(LETTERS, vapply(LETTERS, function(x) paste0(x, LETTERS), character(26)))


# Some mathematical functions
#' @import stats
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
#' @import stats
erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)


#' Determine individual fuzzifier values
#' 
#' This function calculated the values of the fuzzifier from a) the dimensions 
#' of the considered data set and b) 
#' from the individual feature standard deviations.
#' 
#' @param dims vector of two integers containing the dimensions of the data matrix for the clustering
#' @param NClust Number of cluster for running vsclust on (does no influence the calculation of `mm`)
#' @param Sds individual standard deviations, set to 1 if not available
#' @return list of `m`: individual fuzzifiers, `mm`: standard fuzzifier for fcm clustering when not using vsclust algorithm
#' @examples
#' # Generate some random data
#' data <- matrix(rnorm(1:1000), nrow=100)
#' # Estimate fuzzifiers
#' fuzz_out <- determine_fuzz(dim(data), 1)
#' # Run clustering
#' clres <- vsclust_algorithm(data, centers=10, m=fuzz_out$mm)
#' @export
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
determine_fuzz <- function(dims, NClust, Sds = 1) {
  D<-dims[2]
  d<-sqrt(D/2)
  mm <- 1 + (1418/dims[1] + 22.05)* dims[2]^(-2) + (12.33/dims[1] + 0.243) * 
    dims[2]^(-0.0406*log(dims[1])-0.1134)
  
  ### d_i and d_t
  difunc <- function(c,D) { x <- 0:c; sum(choose(c,x)/(x*D+1)*(-1)^x) }
  
  di <- difunc(NClust,D)  /sqrt(pi) * gamma(D/2+1)^(1/D)
  dt <- (NClust)^(-1/D)
  p <- dnorm(di,0,Sds) * (1-dnorm(dt,0,Sds))^(NClust-1)
  
  m <- mm + p*mm*(D/3-1)
  m[m==Inf]<-0
  m[m==0]<-NA
  m[is.na(m)]<-mm*10
  
  ## If m for highest Sd is mm then all = mm
  if (m[which.max(Sds)]== mm) 
    m[seq_len(length(m))] <- mm
  
  return(list (m=m, mm=mm))
  
}


#' Determine optimal cluster number from validity index
#' 
#' Calculated the optimal number from expected behavior of the indices. This would be a large decay for
#' the Minimum Centroid Distance and a minimum for the Xie Beni index
#' 
#' @param ClustInd Output from estimClustNum providing the calculated cluster validity indices
#' @param index Either "MinCentroidDist" or "XieBeni"
#' @param method Either "VSClust" or "FCM" for standard fuzzy c-means clustering
#' @return optimal cluster number
#' @examples
#'   data("artificial_clusters")
#'   dat <- averageCond(artificial_clusters, 5, 10)
#'   dat <- scale(dat)
#' dat <- cbind(dat, 1)
#' ClustInd <- estimClustNum(dat, 10)
#' optimalClustNum
#' @export
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
optimalClustNum <- function(ClustInd, index="MinCentroidDist", method="VSClust") {
  allowedInd <- c("XieBeni", "MinCentroidDist")
  allowedMethod <- c("FCM", "VSClust")
  if (!any(index == allowedInd)) {
    stop(paste("index needs to be one of", paste(allowedInd, collapse=" ")))
  }
  if (!any(method == allowedMethod)) {
    stop(paste("method needs to be one of", paste(allowedMethod, collapse = " ")))
  }
  
  tClustInd <- ClustInd[, grep(index, colnames(ClustInd))]
  tClustInd <- tClustInd[, grep(method, colnames(tClustInd))]
  opt_val <- NULL
  if (length(tClustInd) < 3) 
    stop("Minimal length of ClustInd vector is 3")
  if (index == "MinCentroidDist") {
    opt_val <- which.max(tClustInd[seq_len(length(tClustInd)-1)]-tClustInd[2:length(tClustInd)])
  } else if (index == "XieBeni") {
    opt_val <- which.min(tClustInd[seq_len(length(tClustInd))])
    
  }
  return(as.numeric(opt_val + 2))
}


#' Xie Beni Index of clustering object
#' 
#' Calculate the Xie Beni index for validity of the cluster number in clustering results from running fuzzy c-means or vsclust
#' original publication: 
#' 
#' @references Xie X.L., Beni G. (1991). A validity measure for fuzzy clustering, IEEE Transactions on Pattern Analysis and Machine Intelligence, 13, 841-847.
#' @param clres Output from clustering. Either fclust object or list containing the objects for `membership` and cluster `centers`
#' @param m Fuzzifier value
#' @return Xie Beni index
#' @examples
#' # Generate some random data
#' data <- matrix(rnorm(1:1000), nrow=100)
#' # Run clustering
#' clres <- vsclust_algorithm(data, centers=10, m=1.5)
#' # Calculate Xie-Beni index from results
#' cvalidate.xiebeni(clres, 1.5)
#' @export
cvalidate.xiebeni <- function(clres,m) {                                                
  xrows <- dim(clres$me)[1]                                                
  minimum <- -1                                                            
  error <- clres$within                                                    
  ncenters <- dim(clres$centers)[1]                                        
  for (i in seq_len(ncenters - 1)) {                                            
    for (j in (i + 1):ncenters) {                                        
      diff <- clres$ce[i, ] - clres$ce[j, ]                            
      diffdist <- t(diff) %*% t(t(diff))                               
      if (minimum == -1)                                               
        minimum <- diffdist                                            
      if (diffdist < minimum)                                          
        minimum <- diffdist                                            
    }                                                                    
  }                                                                        
  xiebeni <- error/(xrows * minimum)                                       
  return(xiebeni)                                                          
}

#' Run the vsclust clustering algorithm
#' 
#' This function calls the c++ implementation of the vsclust algorithm, being an extension of fuzzy c-means
#' clustering with additional variance control and capability to run on data with missing values
#'
#' @param x a numeric data matrix 
#' @param centers Either numeric for number of clusters or numeric matrix with center coordinates
#' @param iter.max Numeric for maximum number of iterations
#' @param verbose Verbose information
#' @param dist Distance to use for the calculation. We prefer "euclidean" (default)
#' @param m Fuzzifier value: numeric or vector of length equal to number of rows of x 
#' @param rate.par (experimental) numeric value for punishing missing values
#' @param weights numeric or vector of length equal to number of rows of x 
#' @param control list with arguments to vsclust algorithms (now only cutoff for relative tolerance: reltol)
#' @return list with details about clustering having the objects `centers` (positions of centroids), `size` (feature number per cluster), `cluster` (nearest cluster of each feature), `membership` matrix of membership values, `iter` (number of carried out iterations), `withinerror` (final error from optimization), `call`(call of function)
#' @examples
#' #' # Generate some random data
#' data <- matrix(rnorm(1:1000), nrow=100)
#' # Run clustering
#' clres <- vsclust_algorithm(data, centers=10, m=1.5)
#' head(clres$membership)
#' @export
#' @useDynLib
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
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
    val <- vsclust:::c_plusplus_means(x, centers, weights, m, dist-1, iter.max, 
                                      reltol, verbose, u , 1, iter, NA, rate.par)
    # put modified values in retval
    retval <- list(x=x, xrows = xrows, xcols = xcols, centers = centers,
                   ncenters=ncenters, m = m, dist = dist -1,
                   iter.max = iter.max, reltol = reltol, verbose = verbose, 
                   rate.par = rate.par, u = u, ermin = val, iter = iter)
    
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

#' Plotting vsclust results
#' 
#' This function visualizes the clustered quantitative profiles in multiple figure panels. The parameters allow specifying the main 
#' items like axes labels and color maps. The code is adopted from the MFuzz package.
#' 
#' @param dat a numeric data matrix containing the values used in the clustering
#' @param cl clustering results from vsclust_algorithm or Bestcl object from clustComp function
#' @param mfrow vector of two numbers for the number of rows and colums, figure panels are distributed in the plot
#' @param colo color map to be used (can be missing)
#' @param min.mem filter for showing only features with a higher membership values than this value
#' @param time.labels alternative labels for different conditions
#' @param filename for writing into pdf. Will write on screen when using NA
#' @param xlab Label of x-axis
#' @param ylab Label of y-axis
#' @examples
#' #' # Generate some random data
#' data <- matrix(rnorm(1:5000), nrow=500)
#' # Run clustering
#' clres <- vsclust_algorithm(data, centers=2, m=1.5)
#' mfuzz.plot(data, clres,  mfrow=c(2,3), min.mem=0.0)
#' @return Multiple panels showing expression profiles of clustered features passing the min.mem threshold
#' @export
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.

mfuzz.plot <- function (dat, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels, 
                        filename=NA,xlab="Time",ylab="Expression changes") 
{
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1
  colorindex <- integer(dim(dat)[[1]])
  if (missing(colo)) {
    colo <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", 
              "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", 
              "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", 
              "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", 
              "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", 
              "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", 
              "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", 
              "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", 
              "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", 
              "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", 
              "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", 
              "#FF0060", "#FF0048", "#FF0030", "#FF0018")
  }    else {
    if (colo == "fancy") {
      fancy.blue <- c(c(255:0), rep(0, length(c(255:0))), 
                      rep(0, length(c(255:150))))
      fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
      fancy.red <- c(c(0:255), rep(255, length(c(255:0))), 
                     c(255:150))
      colo <- rgb(blue = fancy.blue/255, green = fancy.green/255, 
                  red = fancy.red/255)
    }
  }
  colorseq <- seq(0, 1, length = length(colo))
  for (j in seq_len(max(clusterindex))) {
    if (sum(clusterindex == j) > 0) {
      tmp <- dat[clusterindex == j, ]
      tmpmem <- memship[clusterindex == j, j]
      if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0) {
        if (!is.na(filename)) {
          pdf(filename, height=3*mfrow[1],width=3*mfrow[2])
        }
        par(mfrow = mfrow, cex=0.5)
        if (sum(clusterindex == j) == 0) {
          ymin <- -1
          ymax <- +1
        }
        else {
          ymin <- min(tmp)
          ymax <- max(tmp)
        }
        plot.default(x = NA, xlim = c(1, dim(dat)[[2]]), 
                     ylim = c(ymin, ymax), xlab = xlab, ylab = ylab, 
                     main = paste("Cluster", j), axes = FALSE)
        if (missing(time.labels)) {
          axis(1, seq_len(dim(dat)[[2]]), c(seq_len(dim(dat)[[2]])))
          axis(2)
        }
        else {
          axis(1, seq_len(dim(dat)[[2]]), time.labels)
          axis(2)
        }
      }
      else {
        if (sum(clusterindex == j) == 0) {
          ymin <- -1
          ymax <- +1
        }
        else {
          ymin <- min(tmp)
          ymax <- max(tmp)
        }
        plot.default(x = NA, xlim = c(1, dim(dat)[[2]]), 
                     ylim = c(ymin, ymax), xlab = xlab, ylab = ylab, 
                     main = paste("Cluster", j), axes = FALSE)
        if (missing(time.labels)) {
          axis(1, seq_len(dim(dat)[[2]]), seq_len(dim(dat)[[2]]))
          axis(2)
        }
        else {
          axis(1, seq_len(dim(dat)[[2]]), time.labels)
          axis(2)
        }
      }
      if (!(sum(clusterindex == j) == 0)) {
        for (jj in seq_len(length(colorseq) - 1)) {
          tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= 
                       colorseq[jj + 1])
          
          if (sum(tmpcol, na.rm=TRUE) > 0) {
            tmpind <- which(tmpcol)
            for (k in seq_len(length(tmpind))) {
              lines(tmp[tmpind[k], ], col = colo[jj])
            }
          }
        }
      }
    }
  }
  if (!is.na(filename))   dev.off()
}

#' Plotting results from estimating the cluster number
#' 
#' This function visualizes the output from estimClustNumber, and there particularly the
#' two validity indices Minimum Centroid Distance and Xie Beni Index.
#' 
#' @param ClustInd Matrix with values from validity indices
#' @return Multiple panels showing expression profiles of clustered features passing the min.mem threshold
#' @examples 
#' data("artificial_clusters")
#' dat <- averageCond(artificial_clusters, 5, 10)
#' dat <- scale(dat)
#' dat <- cbind(dat, 1)
#' ClustInd <- estimClustNum(dat, 10)
#' estimClust.plot(ClustInd)

#' @export
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
estimClust.plot <- function(ClustInd) {
  par(mfrow=c(1,3))
  maxClust <- nrow(ClustInd)+2
  plot(3:maxClust,ClustInd[seq_len(nrow(ClustInd)),"MinCentroidDist_VSClust"],
       col=2 , type="b", 
       main="Min. centroid distance\n(Highest jump is best)",
       xlab="Number of clusters", ylab="Index",
       ylim=c(min(ClustInd[,grep("MinCentroidDist", 
                                 colnames(ClustInd))],na.rm=TRUE),
              max(ClustInd[,grep("MinCentroidDist", 
                                 colnames(ClustInd))],na.rm=TRUE)))
  lines(3:maxClust,ClustInd[seq_len(nrow(ClustInd)),"MinCentroidDist_FCM"],
        col=3,type="b")
  dmindist <- optimalClustNum(ClustInd) 
  points(dmindist,ClustInd[dmindist-2,"MinCentroidDist_VSClust"],pch=15,col=1,cex=2)
  legend("topright",legend = c("VSClust","Standard"), lty=c(1,1),col=2:3)
  grid(NULL,NA,lwd=1,col=1)
  plot(3:maxClust,ClustInd[seq_len(nrow(ClustInd)),"XieBeni_VSClust"], col=2, 
       type="b", main="Xie-Beni index\n(Lowest is best)",
       xlab="Number of clusters", ylab="Index",
       ylim=c(min(ClustInd[,grep("XieBeni", colnames(ClustInd))],na.rm=TRUE),
              max(ClustInd[,grep("XieBeni", colnames(ClustInd))],na.rm=TRUE)))
  lines(3:maxClust,ClustInd[seq_len(nrow(ClustInd)),"XieBeni_FCM"],type="b",col=3)
  dxiebeni <- optimalClustNum(ClustInd, index="XieBeni") 
  points(dxiebeni,ClustInd[dxiebeni-2,"XieBeni_VSClust"],pch=15,col=1,cex=2)
  legend("topright",legend = c("VSClust","Standard"), lty=c(1,1),col=2:3)
  grid(NULL,NA,lwd=1,col=1)
  plot(3:maxClust,ClustInd[seq_len(nrow(ClustInd)),"NumVSClust"], col=2, type="b", 
       main="Total number of assigned features",
       xlab="Number of clusters", ylab="Assigned features",
       ylim=c(min(ClustInd[,grep("Num", colnames(ClustInd))],na.rm=TRUE),
              max(ClustInd[,grep("Num", colnames(ClustInd))],na.rm=TRUE)))
       lines(3:maxClust,ClustInd[seq_len(nrow(ClustInd)),"NumFCM"],type="b",col=3)
       legend("topright",legend = c("VSClust","Standard"), lty=c(1,1),col=2:3)
       # finally plot
       p <- recordPlot()
}

#' arrange cluster member numbers from largest to smallest
#' @param Bestcl fclust object
#' @param NClust Number of clusters
#' @importFrom matrixStats rowMaxs
SwitchOrder <- function(Bestcl,NClust) {
  switching <- as.numeric(names(sort(table(Bestcl$cluster[rowMaxs(Bestcl$membership)>0.5]),
                                     decreasing=TRUE)))
  if (length(switching)<NClust) 
    switching <- c(switching,which(!((seq_len(NClust)) %in% switching)))
  switching2 <- seq_len(NClust)
  names(switching2) <- switching
  tBest <- Bestcl
  tBest$centers <- Bestcl$centers[switching,]
  rownames(tBest$centers) <- seq_len(NClust)
  tBest$size <- Bestcl$size[switching]
  tBest$cluster <- switching2[as.character(Bestcl$cluster)]
  names(tBest$cluster) <- names(Bestcl$cluster)
  tBest$membership <- Bestcl$membership[,switching]
  colnames(tBest$membership) <- seq_len(NClust)
  tBest
}

#' Function to run clustering with automatic fuzzifier settings (might become obsolete)
#' 
#' Run original fuzzy c-means and vsclust for a number of clusters and the given data set including data pre-processing and
#' automatic setting of the data-dependent parameters like the lower limit of the fuzzifier.
#' 
#' @param dat a numeric data matrix 
#' @param NSs number of clusterings runs with different random seeds
#' @param NClust Number of clusters
#' @param Sds Standard deviation of features (either vector of the same length as features numbers in matrix or single value)
#' @param cl object of class `cluster` or `SOCKcluster` to specify environment for parallelization
#' @return List containing the objects 
#' @return `indices` containing minimum centroid distance and Xie-Beni index for both clustering methods
#' @return `Bestcl` optimal vsclust results (variance-sensitive fcm clustering)
#' @return `Bestcl2` optimal fuzzy c-means restults 
#' @return `m` vector of individual fuzzifer values per feature
#' @return `withinerror` final optimization score for vsclust
#' @return `withinerror2` final optimization score for fuzzy c-means clustering
#' @examples
#' #' # Generate some random data
#' data <- matrix(rnorm(1:1000), nrow=100)
#' # Run clustering
#' cl <- parallel::makePSOCKcluster(1, nnodes=1)
#' ClustCompOut <- ClustComp(data, cl=cl, NClust=10, Sds=1)
#' barplot(ClustCompOut$indices)
#' @import parallel
#' @import stats
#' @export
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
ClustComp <- function(dat,NSs=10,NClust=NClust,Sds=Sds, cl=parallel::makePSOCKcluster(1)) {
  fuzz_out <- determine_fuzz(dim(dat), NClust, Sds)
  m <- fuzz_out$m
  mm <- fuzz_out$mm
  
  clusterExport(cl=cl,varlist=c("dat","NClust","m","mm"),envir=environment())
  cls <- parLapply(cl,seq_len(NSs), function(x) vsclust_algorithm(dat,NClust,
                                                                  m=m,verbose=F,
                                                                  iter.max=1000))
  # cls <- lapply(1:NSs, function(x) vsclust_algorithm(tData,NClust,m=m,verbose=F,iter.max=1000))  #print(cls[[1]])
  Bestcl <- cls[[which.min(lapply(cls,function(x) x$withinerror))]]
  cls <- parLapply(cl,seq_len(NSs), function(x) vsclust_algorithm(dat,NClust,
                                                                  m=mm,verbose=F,
                                                                  iter.max=1000))
  Bestcl2 <- cls[[which.min(lapply(cls,function(x) x$withinerror))]]
  # stopCluster(cl)
  
  # return validation indices
  list(indices=c(min(dist(Bestcl$centers)),cvalidate.xiebeni(Bestcl,mm),
                 min(dist(Bestcl2$centers)),cvalidate.xiebeni(Bestcl2,mm)),
       Bestcl=Bestcl,Bestcl2=Bestcl2,m=m,withinerror=Bestcl$withinerror,
       withinerror2=Bestcl2$withinerror) 
}

#' Paired statistical testing
#' 
#' Statistical testing and variance estimation in multi-dimensional data set. given by a matrix. This functions runs LIMMA paired tests and 
#' calculated the shrunken variance estimates.
#' 
#' @param Data a numeric data matrix with columns as samples. Different experimental conditions are grouped together in their replicates. The number of samples per group needs to be identical
#' @param NumCond Number of different experimental conditions
#' @param NumReps Number of replicates per experimental condition
#' @return List containing the objects 
#' @return `qvalues` false discovery rates after correction for multiple testing (`qvalue` method from `qvalue` library)
#' @return `Sds` General standard deviation within replicates after using shrinkage (eBayes) by LIMMA
#' @examples
#' #' # Generate some random data with three different experimental conditions
#' data <- matrix(rnorm(1:1500), nrow=100)
#' # Run statistical testing
#' stat_out <- SignAnalysisPaired(data, 3, 5)
#' # Histogram of qvalues comparing the second to the first condition
#' hist(stat_out$qvalues[,1], 50, xlab="q-values")
#' @import limma
#' @import stats
#' @importFrom qvalue qvalue
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
#' @export
SignAnalysisPaired <- function(Data,NumCond,NumReps) {
  ##########################################################
  # significance analysis
  MAData<-Data[,2:(NumCond)]-Data[,1]
  for (i in seq_len(NumReps-1))
    MAData<-cbind(MAData,Data[,(i*NumCond+1)+seq_len(NumCond-1)]-Data[,(i*NumCond+1)])
  rownames(MAData)<-rownames(Data)
  MAReps<-rep(seq_len(NumCond-1),NumReps)
  if (is.null(rownames(MAData))) rownames(MAData) <- paste0("feature",seq_len(nrow(MAData)))
  ##limma with ratios
  design<-plvalues<-NULL
  for (c in (seq_len(NumCond-1))) {
    design<-cbind(design,as.numeric(MAReps==c))
  }
  lm.fittedMA <- lmFit(MAData,design)
  lm.bayesMA<-eBayes(lm.fittedMA)
  topTable(lm.bayesMA)
  plvalues <- lm.bayesMA$p.value
  qvalues <- matrix(NA,nrow=nrow(plvalues),ncol=ncol(plvalues),dimnames=dimnames(plvalues))
  # qvalue correction
  for (i in seq_len(ncol(plvalues))) {
    tqs <- tryCatch(qvalue(na.omit(plvalues[,i]))$qvalues, 
                    error = function(e) NULL)
    if (length(tqs) >0) {
      qvalues[names(tqs),i] <- tqs
    }
    else {
      qvalues[names(tqs),i] <- NA
    }
  }
  
  return(list(qvalues=qvalues,Sds=sqrt(lm.bayesMA$s2.post)))
}

#' Unpaired statistical testing
#' 
#' Statistical testing and variance estimation in multi-dimensional data set. given by a matrix. This functions runs LIMMA paired tests and 
#' calculated the shrunken variance estimates.
#'
#' @param Data a numeric data matrix with columns as samples. Different experimental conditions are grouped together in their replicates. The number of samples per group needs to be identical
#' @param NumCond Number of different experimental conditions
#' @param NumReps Number of replicates per experimental condition
#' @return List containing the objects 
#' @return `pvalues` p-values before correction for multiple testing
#' @return `qvalues` false discovery rates after correction for multiple testing (`qvalue` method from `qvalue` library)
#' @return `Sds` General standard deviation within replicates after using shrinkage by LIMMA
#' @examples
#' #' # Generate some random data
#' data <- matrix(rnorm(1:1000), nrow=100)
#' # Run statistical testing
#' stat_out <- SignAnalysis(data, 2, 5)
#' # Histogram of qvalues (no significant events)
#' hist(stat_out$qvalues, 50, xlab="q-values")
#' @import limma
#' @import stats
#' @importFrom qvalue qvalue
#' @export
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
SignAnalysis <- function(Data,NumCond,NumReps) {
  ##########################################################
  if (is.null(rownames(Data))) rownames(Data) <- paste0("feature",seq_len(nrow(Data)))
  # significance analysis
  Reps <- rep(seq_len(NumCond),NumReps)
  design <- model.matrix(~0+factor(Reps-1))
  colnames(design)<-paste("i",c(seq_len(NumCond)),sep="")
  contrasts<-NULL
  First <- 1
  for (i in (seq_len(NumCond))[-First]) contrasts<-append(contrasts,
                                                          paste(colnames(design)[i],"-",colnames(design)[First],sep=""))
  contrast.matrix<-makeContrasts(contrasts=contrasts,levels=design)
  lm.fitted <- lmFit(Data,design)
  
  lm.contr <- contrasts.fit(lm.fitted,contrast.matrix)
  lm.bayes<-eBayes(lm.contr)
  #topTable(lm.bayes)
  plvalues <- lm.bayes$p.value
  qvalues <- matrix(NA,nrow=nrow(plvalues),ncol=ncol(plvalues),
                    dimnames=dimnames(plvalues))
  # qvalue correction
  for (i in seq_len(ncol(plvalues))) {
    tqs <- tryCatch(qvalue(na.omit(plvalues[,i]))$qvalues, 
                    error = function(e) NULL)
    # print(tqs)
    if (length(tqs) > 0) {
      qvalues[names(tqs),i] <- tqs
    }
    else {
      qvalues[names(tqs),i] <- NA
    }
  }
  return(list(pvalues=plvalues,qvalues=qvalues,Sds=sqrt(lm.bayes$s2.post)))
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
#' the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. 
#' doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters 
#' for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. 
#' doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.

calcBHI <- function(Accs,gos) {
  ## enrichment might not yield all GO terms! This could lead to problems
  BHI <- sumcomb <- vector("integer",length(Accs))
  names(BHI) <- names(Accs)
  goData <- as.data.frame(gos@compareClusterResult)
  for (i in names(gos@geneClusters)) {
    genes <- Accs[[i]]
    ispair <- matrix(F,length(genes),length(genes),dimnames = list(rows=genes,cols=genes))
    sumcomb[i] <- choose(length(genes),2)
    clgroup <- goData[goData$Cluster==i,"geneID"]
    for (j in seq_len(length(clgroup))) {
      tgenes <- genes[na.omit(match(unlist(strsplit(clgroup[j],"/"),use.names=F),as.character(genes)))]
      ltgenes <- length(tgenes)
      if (ltgenes > 1) {
        for (i1 in tgenes[seq_len(ltgenes-1)]) {
          ttgene <- tgenes[(which(i1==tgenes)+1):ltgenes]
          ispair[i1,ttgene] <- ispair[ttgene,i1] <- TRUE
        }
        
      }
      # print(tgenes)
    }
    BHI[i] <- BHI[i] + sum(ispair)/2
  }
  # Rprof()
  # print(summaryRprof(tmp))
  sum(BHI)/sum(sumcomb)
}

#' Calculate mean over replicates
#' 
#' Simple method to calculate the means for each feature across its replicates
#' 
#' @param data Matrix of data frame with numerical values. Columns corresponds to samples
#' @param NumReps Number of replicates per experimental condition
#' @param NumCond Number of different experimental conditions
#' @return Matrix of data frame with averaged values over replicates for each conditions
#' @examples
#' data <- matrix(rnorm(1000), nrow=100)
#' av_data <- averageCond(data, NumCond=2, NumReps=5)
#' @export
averageCond <- function(data, NumReps, NumCond) {
  # Calculates means over replicates
  tdat<-rowMeans(data[,seq(1,NumReps*NumCond,NumCond)],na.rm=TRUE)
  for (i in 2:NumCond) {
    tdat<-cbind(tdat,rowMeans(data[,seq(i,NumReps*NumCond,NumCond)],na.rm=TRUE))
  }
  colnames(tdat)<-paste("Mean of log ",LETTERS702[seq_len(NumCond)],sep="")
  tdat  
  
}

#' Visualize using principal component analysis (both loadings and scoring) including the variance from the replicates
#' 
#' The loading plot shows all features and their scaled variance. This provides an idea of the intrinsic noise in the data.
#' 
#' @param data Matrix of data frame with numerical values. Columns corresponds to samples
#' @param NumReps Number of replicates per experimental condition
#' @param NumCond Number of different experimental conditions
#' @param Sds Standard deviation for each features. Usually using the one from LIMMA
#' @return Loading and scoring plots that include feature variance
#' @examples
#' data <- matrix(rnorm(1000), nrow=100)
#' pcaWithVar(data, NumCond=2, NumReps=5, Sds=1)
#' @import stats
#' @import graphics
#' @importFrom shiny need
#' @importFrom grDevices rainbow
#' @export
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
pcaWithVar <- function(data, NumReps, NumCond, Sds=1) {
  # Remove columns with only 5% of the data
  plab <- rep(seq_len(NumCond),NumReps)
  plab <- plab[colSums(!is.na(data))>0.05*nrow(data)]
  pcaDat <- data
  if (ncol(pcaDat) != NumCond * NumReps)
    stop("Wrong number of conditions and/or replicates!")
  pcaDat <- data[,colSums(!is.na(data))>0.05*nrow(data)]
  pcaDat <- (pcaDat[complete.cases(pcaDat),])
  validate(need(length(pcaDat)> 0, "Principal component analysis not shown as too many missing values"))      
  validate(need(nrow(pcaDat)> 10, "Principal component analysis not shown as too many missing values"))      
  pca<-prcomp(pcaDat,scale=TRUE,retx=TRUE)
  # scores <- pca$x
  # loadings <- pca$rotation
  scores <- pca$rotation
  loadings <- pca$x
  par(mfrow=c(1,2))
  # Set missing values to maximum std. dev.
  Sds[is.na(Sds)] <- max(Sds, na.rm=TRUE)
  # Scaling suitable for visualization
  Sds <- sqrt(Sds)
  plot(loadings,cex=Sds,pch=16,
       col=paste("#000000",sprintf("%02X",as.integer(255/max(1/Sds)/Sds)),sep=""))
  title(main="Principal component analysis of data set (loadings)",sub="The point size corresponds to the estimated standard deviation")
  plot(scores,pch=19,col=rainbow(NumCond)[plab])
  title(main="Principal component analysis of data set (scores)",sub="Colors denote different conditions")
  legend("topright",paste("Condition",seq_len(NumCond)),pch=rep(19,NumCond),
         col=rainbow(NumCond)[seq_len(NumCond)])
}

### Wrapper functions

#' Wrapper for statistical analysis
#' 
#' Prepare data for running vsclust clustering. 
#' This includes visualization runnig the functions for the principal component analysis and its visualization, 
#' statistical testing with LIMMA, as well as scaling and filtering of missing values
#' @param dat matrix or data frame of numerical data. Columns are samples. Replicates are grouped (i.e. A1, B1, C1, A2, B2, C2) when letters denote conditions and numbers the replicates. In case of `isStat=F`, you need a last column for the standard deviations
#' @param NumReps Number replicates in the data
#' @param NumCond Number of different experimental conditions. The total number of columns needs to be NumReps*NumCond 
#' @param isPaired Boolean for running paired or unpaired statistical tests
#' @param isStat Boolean for whether to run statistical test or each column corresponds to a different experimental conditions. Then this function reads feature standard deviations from data frame from the last column
#' @return list with the items `dat` (data matrix of features averaged over replicates and last column with their standard deviations), `qvals` FDRs from the statistical tests (each conditions versus the first), `StatFileOut` all of before for saving in file
#' @examples
#' data <- matrix(rnorm(2000), nrow=200)
#' stats <- PrepareForVSClust(data, 5, 2, isStat=TRUE)
#' 
#' @import stats
#' @importFrom matrixStats rowSds
#' @importFrom shiny validate
#' @export
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
PrepareForVSClust <- function(dat, NumReps, NumCond, isPaired=F, isStat) {
  qvals <- statFileOut <- Sds <- NULL
  tdat <- NULL
  
  # Run statistical testing
  if (isStat) {
    if(ncol(dat)!=NumReps*NumCond)
      stop("Number of data columns must correspond to product of conditions and replicates!")
    if (isPaired) {
      ttt <- SignAnalysisPaired(dat,NumCond, NumReps)
    } else {
      ttt <- SignAnalysis(dat,NumCond, NumReps)
    }
    
    Sds <- ttt$Sds
    qvals <- ttt$qvalues
    colnames(qvals) <- paste("qvalue ",LETTERS702[2:(NumCond)],"vsA",sep="")
    
    tdat <- averageCond(dat, NumReps, NumCond)
    
  } else {
    Sds <- dat[,ncol(dat)]
    tdat <- dat[,seq_len(ncol(dat)-1)]
    NumReps <- 1
    NumCond <- ncol(dat)-1
    dat <- tdat
  }
  
  if (isStat) {
    statFileOut <- cbind(tdat,Sds,qvals)
  } else {
    statFileOut <- cbind(tdat,Sds)
  }
  
  pcaWithVar(dat, NumReps, NumCond, Sds / rowSds(tdat, na.rm=TRUE))
  
  ## Preparing output
  Out <- list(dat=cbind(tdat,Sds), qvals=qvals, statFileOut=statFileOut)
  Out
  
}

#' Wrapper for estimation of cluster number
#' 
#' This runs the clustering for different numbers of clusters, and estimates the most suitable numbers from applying
#' the minimum centroid distance and the Xie Beni index. Multi-threading is used to shorten the computation times. 
#' Given the hierarchical structure of many data sets, the resulting 
#' numbers are suggestions. Inspection of the here plotted indices help to determine alternative cluster numbers, 
#' given by a strong decay of the minimum centroid distance and/or a low value of the Xie Beni index.
#' 
#' @param dat matrix of features averaged over replicates. The last column contains their standard deviation
#' @param maxClust Maximal number of cluster. The minimum is 3
#' @param cores The number of threads to be used for parallelisation
#' @return list with the items `ClustInd`: list of clustering objects for each number of clusters, `p` plot object with plots for validity indices, `numclust` optimal cluster number according to "minimum centroid distance"
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
estimClustNum<- function(dat, maxClust=25, cores=1) {
  ClustInd<-matrix(NA,nrow=maxClust-2,ncol=6)
  if (is.null(rownames(dat))) rownames(dat) <- seq_len(nrow(dat))
  tData <- dat[,seq_len(ncol(dat)-1)]
  colnames(tData)<-NULL
  
  # define parallelization
  cl <- makeCluster(cores)
  clusterExport(cl=cl,varlist=c("vsclust_algorithm"),envir=environment())
  clusterEvalQ(cl=cl, library(vsclust))
  
  # Standardise
  sds <- dat[rownames(tData), ncol(dat)]
  # scale standard deviations by the ones in the actual data to cope for the following standardization 
  sds <- sds / (rowSds(as.matrix(tData),na.rm=TRUE))
  tData <- t(scale(t(tData)))
  
  multiOut <- lapply(3:maxClust,function(x) {    
    if (!is.null(getDefaultReactiveDomain())) {
      incProgress(1, detail = paste("Running cluster number",x))
    } else {
      message(paste("Running cluster number",x))
    }
    clustout <- ClustComp(tData,NClust=x,Sds=sds,NSs=16, cl=cl)
    c(clustout$indices,sum(rowMaxs(clustout$Bestcl$membership)>0.5),
      sum(rowMaxs(clustout$Bestcl2$membership)>0.5))
  })
  
  stopCluster(cl)
  
  for (NClust in 3:maxClust) 
    ClustInd[NClust-2,] <- multiOut[[NClust-2]]
  rownames(ClustInd) <- paste0("num_clust_",3:maxClust)
  colnames(ClustInd) <- c("MinCentroidDist_VSClust", "XieBeni_VSClust", "MinCentroidDist_FCM", "XieBeni_FCM", "NumVSClust", "NumFCM")
  
  # Output
  ClustInd
}

#' Wrapper for running cluster analysis
#' 
#' This function runs the clustering and visualizes the results.
#' 
#' @param dat matrix or data frame with feature values for different conditions
#' @param NClust Number of cluster for running the clustering
#' @param proteins vector with additional feature information (default is NULL) to be added to the results
#' @param VSClust boolean. TRUE for running the variance-sensitive clustering. Otherwise, the function will call standard fuzzy c-means clustering
#' @param cores Number of threads for the parallelization
#' @return list with the items `dat`(the original data), `Bestcl` clustering results (same as from vsclust_algorithm), `p` (plot object with mfuzz plots), `outFileClust`(suitable matrix with complete information) , `ClustInd` (information about being member of any cluster, feature needs on membership values > 0.5)
#' @examples
#' data(iris)
#' data <- cbind(iris[,1:4],1)
#' clust_out <- runClustWrapper(data, NClust=3, cores=1)
#' clust_out$p
#' @import parallel
#' @import graphics
#' @importFrom grDevices recordPlot
#' @importFrom shiny getDefaultReactiveDomain incProgress
#' @importFrom matrixStats rowMaxs
#' @export
runClustWrapper <- function(dat, NClust, proteins=NULL, VSClust=TRUE, cores) {
  tData <- dat[,seq_len(ncol(dat)-1)]
  sds <- dat[,ncol(dat)]
  
  #Standardize
  # scale standard deviations by the ones in the actual data to cope for the following standardization 
  sds <- sds / rowSds(as.matrix(tData),na.rm=TRUE)
  tData <- t(scale(t(tData)))
  if (is.null(rownames(tData))) {
    rownames(tData) <- seq_len(nrow(tData))
  }
  cl <- makeCluster(cores)
  clusterExport(cl=cl,varlist=c("vsclust_algorithm"),envir=environment())
  clusterEvalQ(cl=cl, library(vsclust))
  
  clustout <- ClustComp(tData,NClust=NClust,Sds=sds, NSs=16, cl=cl)
  stopCluster(cl)
  
  
  if (VSClust) {
    Bestcl <- clustout$Bestcl
  } else {
    Bestcl <- clustout$Bestcl2
  }
  Bestcl <- SwitchOrder(Bestcl,NClust)
  
  # sorting for membership values (globally)
  Bestcl$cluster <- Bestcl$cluster[order(rowMaxs(Bestcl$membership,na.rm=TRUE))]
  Bestcl$membership <- Bestcl$membership[order(rowMaxs(Bestcl$membership,na.rm=TRUE)),]
  tData <- tData[names(Bestcl$cluster),]
  
  if (!is.null(getDefaultReactiveDomain()))
    incProgress(0.7, detail = paste("Plotting",NClust))
  
  # graphics.off() ## clean up device
  par(lwd=0.25)
  oldmar <- par("mar")
  par(mar=c(2,2,3,3),mgp=c(2,1,0))
  par(mar=par("mar")/max(1,NClust/20))
  
  mfuzz.plot(tData,cl=Bestcl,mfrow=c(round(sqrt(NClust)),ceiling(sqrt(NClust))),
             min.mem=0.5,colo="fancy")
  p <- recordPlot()
  # par(lwd=1,mar=oldmar)
  
  colnames(Bestcl$membership) <- paste("membership of cluster",colnames(Bestcl$membership))
  outFileClust <- tData
  if (!is.null(proteins)) {
    outFileClust <- cbind(outFileClust,names=as.character(proteins[rownames(outFileClust)]))
  }
  
  rownames(Bestcl$centers) <- paste("Cluster",rownames(Bestcl$centers))
  ClustInd <- as.data.frame(table(Bestcl$cluster[rowMaxs(Bestcl$membership)>0.5]))
  if (ncol(ClustInd) == 2)
    colnames(ClustInd) <- c("Cluster","Members")
  else 
    ClustInd <- cbind(seq_len(max(Bestcl$cluster)),rep(0,max(Bestcl$cluster)))
  
  ## Output
  Out <- list(dat=tData, Bestcl=Bestcl, p=p, outFileClust=outFileClust, ClustInd=ClustInd)
  return(Out)
}

#' Wrapper for functional enrichment
#' 
#' The functional analysis uses the libarary RDAVIDWebService and thus might become obsolete as that library is not supported anymore
#' The user can select different ID types and different enrichment categories like GO terms and pathways. 
#' Allowed ID types: 
#' "AFFYMETRIX_3PRIME_IVT_ID", 
#' "AFFYMETRIX_EXON_GENE_ID", "AGILENT_CHIP_ID", 
#' "AGILENT_ID", "AGILENT_OLIGO_ID", "APHIDBASE_ID", "BEEBASE_ID", 
#' "BEETLEBASE_ID", "BGD_ID", "CGNC_ID", "CRYPTODB_ID", "DICTYBASE_ID", "ENSEMBL_GENE_ID", 
#' "ENSEMBL_TRANSCRIPT_ID", "ENTREZ_GENE_ID", "GENOMIC_GI_ACCESSION", "FLYBASE_GENE_ID", "GENBANK_ACCESSION",
#' "GENPEPT_ACCESSION", "LOCUS_TAG", "ILLUMINA_ID", "MGI_ID", "MIRBASE_ID",
#' "OFFICIAL_GENE_SYMBOL", "PFAM_ID", "PIR_ID", "PROTEIN_GI_ACCESSION", "MRNA_GI_ACCESSION",
#' "REFSEQ_GENOMIC", "REFSEQ_MRNA", "REFSEQ_PROTEIN", "REFSEQ_RNA", 
#' "RGD_ID", "SGD_ID", "TAIR_ID", "UCSC_GENE_ID", "UNIGENE", 
#' "UNIPROT_ACCESSION", "UNIPROT_ID", "UNIREF100_ID", "WORMBASE_GENE_ID", 
#' "WORMPEP_ID", "ZFIN_ID"
#' Allowed enrichment categories:
#' "GOTERM_MF_ALL","GOTERM_BP_ALL",
#' "GOTERM_CC_ALL","GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"),
#' "KEGG"="KEGG_PATHWAY","PANTHER_PATHWAY","REACTOME_PATHWAY","BBID","BIOCARTA",
#' "DIP","MINT","INTACT","BIOGRID_INTERACTION","GAD_DISEASE","GAD_DISEASE_CLASS",
#' "OMIM_DISEASE","INTERPRO","PROSITE","PFAM","SMART","PRODOM","PIR_SUPERFAMILY"
#' 
#' @param cl clustering results (either directly from vsclust_algorithm or as `Bestcl` object from ClustComp or runClustWrapper)
#' @param protnames vector providing the corresponding gene/protein names of the features (set to NULL for directly using the feature names (default))
#' @param idtypes type of IDs for features given by genes/proteins (generic gene names are not working)
#' @param infosource Type of gene annotation (e.g. KEGG_PATHWAY)
#' @return plot object to be able to pass the figures to e.g. shiny
#' @examples
#' \donttest{ library(clusterProfiler)
#'  data(gcSample)
#' data <- cbind(matrix(rnorm(2000), nrow=500), sds=1)
#' # Adding an artificial cluster
#' data[,c(1,3)] <- data[,c(1,3)] + rep(c(1,-1),250)
#' rownames(data) <- gcSample[[7]][1:500]
#'
#' clust_out <- runClustWrapper(data, NClust=2, cores=1)
#' # Taking some gene names from example data set
#' # This function calls the DAVID web service and thus will require an internet connection
#' enrich_out <- runFuncEnrich(clust_out$Bestcl, NULL, "ENTREZ_GENE_ID", "GOTERM_MF_ALL")
#' dotplot(enrich_out$fullFuncs)
#' }
#' @import graphics
#' @importFrom clusterProfiler compareCluster
#' @importFrom shiny need
#' @importFrom matrixStats rowMaxs
#' @export
runFuncEnrich <- function(cl, protnames=NULL, idtypes, infosource) {
  Accs <- list()
  for (c in seq_len(max(cl$cluster))) {
    cname <- paste("Cluster",c, sep="_")
    Accs[[cname]] <- names(which(cl$cluster==c & rowMaxs(cl$membership)>0.5))
    
    Accs[[cname]] <- Accs[[cname]][Accs[[cname]]!=""]
    if (length(Accs[[cname]])>0) {
      if (!is.null(protnames)) {
        Accs[[cname]] <- as.character(protnames[Accs[[cname]]])
      }
      
      Accs[[cname]] <- sub("-[0-9]","",Accs[[cname]])
    }
  }
  # TODO? add extraction of multiple accession numbers
  Accs <- lapply(Accs,function(x) unique(ifelse(is.na(x),"B3",x)))
  Accs <- Accs[lapply(Accs,length)>0]
  x <- NULL
  try(x <- compareCluster(Accs, fun="enrichDAVID", annotation=infosource,
                          idType=idtypes,
                          david.user = "veits@bmb.sdu.dk"))
  validate(need(!is.null(x),"No result. Wrong ID type?"))
  if (!is.null(getDefaultReactiveDomain()))
    incProgress(0.7, detail = "received")
  message("got data from DAVID\n")
  x@compareClusterResult <- cbind(x@compareClusterResult,
                                  log10padval=log10(x@compareClusterResult$p.adjust))
  y <- new("compareClusterResult",compareClusterResult=x@compareClusterResult)
  if (length(unique(y@compareClusterResult$ID)) > 20) {
    message("Reducing number of DAVID results\n")
    y@compareClusterResult <- y@compareClusterResult[
      order(y@compareClusterResult$p.adjust)[seq_len(20)],]
    
    y@compareClusterResult$Cluster <- as.character(y@compareClusterResult$Cluster)
  }
  
  BHI <- calcBHI(Accs,x)
  return(list(fullFuncs=x, redFuncs=y, BHI=BHI))
  
}

#' Run VSClust as Shiny app
#' 
#' You will get the full functionality of the VSClust workflow with multiple visualizations and downloads
#' 
#' @return The shiny app should open in a browser or in RStudio.
#' @examples
#' runVSClustApp()
#' @export
#' @references 
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#' 
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#' 
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
runVSClustApp <- function() {
  shiny::runApp(system.file("shiny/", package="vsclust"))
}