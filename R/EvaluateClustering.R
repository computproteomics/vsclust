#' 
#' Xie Beni Index of clustering object
#'
#' Calculate the Xie Beni index for validity of the cluster number in clustering 
#' results from running fuzzy c-means or vsclust
#' original publication:
#'
#' @references Xie X.L., Beni G. (1991). A validity measure for fuzzy 
#' clustering, IEEE Transactions on Pattern Analysis and Machine Intelligence, 
#' 13, 841-847.
#' @param clres Output from clustering. Either fclust object or list containing 
#' the objects for `membership` and cluster `centers`
#' @param m Fuzzifier value
#' @return Xie Beni index
#' @examples
#' # Generate some random data
#' data <- matrix(rnorm(seq_len(1000)), nrow=100)
#' # Run clustering
#' clres <- vsclust_algorithm(data, centers=5, m=1.5)
#' # Calculate Xie-Beni index from results
#' cvalidate.xiebeni(clres, 1.5)
#' @export
cvalidate.xiebeni <-
  function(clres, m) {
    xrows <-
      dim(clres$me)[1]
    minimum <-
      -1
    error <-
      clres$within
    ncenters <-
      dim(clres$centers)[1]
    for (i in seq_len(ncenters - 1)) {
      for (j in seq(i + 1,ncenters, 1)) {
        diff <- clres$ce[i,] - clres$ce[j,]
        diffdist <-
          t(diff) %*% t(t(diff))
        if (minimum == -1)
          minimum <-
          diffdist
        if (diffdist < minimum)
          minimum <-
          diffdist
      }
    }
    xiebeni <-
      error / (xrows * minimum)
    return(xiebeni)
  }


#' Plotting vsclust results
#'
#' This function visualizes the clustered quantitative profiles in multiple 
#' figure panels. The parameters allow specifying the main items like axes 
#' labels and color maps. The code is adopted from the MFuzz package.
#'
#' @param dat a numeric data matrix containing the values used in the clustering
#' @param cl clustering results from vsclust_algorithm or Bestcl object from 
#' clustComp function
#' @param mfrow vector of two numbers for the number of rows and colums, figure 
#' panels are distributed in the plot
#' @param colo color map to be used (can be missing)
#' @param minMem filter for showing only features with a higher membership 
#' values than this value
#' @param timeLabels alternative labels for different conditions
#' @param filename for writing into pdf. Will write on screen when using NA
#' @param xlab Label of x-axis
#' @param ylab Label of y-axis
#' @examples
#' #' # Generate some random data
#' data <- matrix(rnorm(seq_len(5000)), nrow=500)
#' # Run clustering
#' clres <- vsclust_algorithm(data, centers=2, m=1.5)
#' mfuzz.plot(data, clres,  mfrow=c(2,3), minMem=0.0)
#' @return Multiple panels showing expression profiles of clustered features 
#' passing the minMem threshold
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
#' for fuzzy c-means cluster analysis. Bioinformatics. 
#' 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 
#' 29. PMID: 20880957.

mfuzz.plot <-
  function (dat,
            cl,
            mfrow = c(1, 1),
            colo,
            minMem = 0,
            timeLabels,
            filename = NA,
            xlab = "Time",
            ylab = "Expression changes")
  {
    clusterindex <- cl[[3]]
    memship <- cl[[4]]
    memship[memship < minMem] <- -1
    colorindex <- integer(dim(dat)[[1]])
    if (missing(colo)) {
      colo <- c(
        "#FF8F00",
        "#FFA700",
        "#FFBF00",
        "#FFD700",
        "#FFEF00",
        "#F7FF00",
        "#DFFF00",
        "#C7FF00",
        "#AFFF00",
        "#97FF00",
        "#80FF00",
        "#68FF00",
        "#50FF00",
        "#38FF00",
        "#20FF00",
        "#08FF00",
        "#00FF10",
        "#00FF28",
        "#00FF40",
        "#00FF58",
        "#00FF70",
        "#00FF87",
        "#00FF9F",
        "#00FFB7",
        "#00FFCF",
        "#00FFE7",
        "#00FFFF",
        "#00E7FF",
        "#00CFFF",
        "#00B7FF",
        "#009FFF",
        "#0087FF",
        "#0070FF",
        "#0058FF",
        "#0040FF",
        "#0028FF",
        "#0010FF",
        "#0800FF",
        "#2000FF",
        "#3800FF",
        "#5000FF",
        "#6800FF",
        "#8000FF",
        "#9700FF",
        "#AF00FF",
        "#C700FF",
        "#DF00FF",
        "#F700FF",
        "#FF00EF",
        "#FF00D7",
        "#FF00BF",
        "#FF00A7",
        "#FF008F",
        "#FF0078",
        "#FF0060",
        "#FF0048",
        "#FF0030",
        "#FF0018"
      )
    }    else {
      if (colo == "fancy") {
        fancyBlue <- c(c(255:0), rep(0, length(c(255:0))),
                        rep(0, length(c(255:150))))
        fanceGreen <-
          c(c(0:255), c(255:0), rep(0, length(c(255:150))))
        fancyRed <- c(c(0:255), rep(255, length(c(255:0))),
                       c(255:150))
        colo <- rgb(blue = fancyBlue / 255,
                    green = fanceGreen / 255,
                    red = fancyRed / 255)
      }
    }
    colorseq <- seq(0, 1, length = length(colo))
    for (j in seq_len(max(clusterindex))) {
      if (sum(clusterindex == j) > 0) {
        # keep matrices in place even if a cluster has just one member
        tmp <- dat[clusterindex == j,, drop = FALSE]
        tmpmem <- memship[clusterindex == j, j]
        if (((j - 1) %% (mfrow[1] * mfrow[2])) == 0) {
          if (!is.na(filename)) {
            pdf(filename,
                height = 3 * mfrow[1],
                width = 3 * mfrow[2])
          }
          par(mfrow = mfrow, cex = 0.5)
          if (sum(clusterindex == j) == 0) {
            ymin <- -1
            ymax <- +1
          }
          else {
            ymin <- min(tmp, na.rm = TRUE)
            ymax <- max(tmp, na.rm = TRUE)
          }
          plot.default(
            x = NA,
            xlim = c(1, dim(dat)[[2]]),
            ylim = c(ymin, ymax),
            xlab = xlab,
            ylab = ylab,
            main = paste("Cluster", j),
            axes = FALSE
          )
          if (missing(timeLabels)) {
            axis(1, seq_len(dim(dat)[[2]]), c(seq_len(dim(dat)[[2]])))
            axis(2)
          }
          else {
            axis(1, seq_len(dim(dat)[[2]]), timeLabels)
            axis(2)
          }
        }
        else {
          if (sum(clusterindex == j) == 0) {
            ymin <- -1
            ymax <- +1
          }
          else {
            ymin <- min(tmp, na.rm = TRUE)
            ymax <- max(tmp, na.rm = TRUE)
          }
          plot.default(
            x = NA,
            xlim = c(1, dim(dat)[[2]]),
            ylim = c(ymin, ymax),
            xlab = xlab,
            ylab = ylab,
            main = paste("Cluster", j),
            axes = FALSE
          )
          if (missing(timeLabels)) {
            axis(1, seq_len(dim(dat)[[2]]), seq_len(dim(dat)[[2]]))
            axis(2)
          }
          else {
            axis(1, seq_len(dim(dat)[[2]]), timeLabels)
            axis(2)
          }
        }
        if (!(sum(clusterindex == j) == 0)) {
          for (jj in seq_len(length(colorseq) - 1)) {
            tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <=
                         colorseq[jj + 1])
            
            if (sum(tmpcol, na.rm = TRUE) > 0) {
              tmpind <- which(tmpcol)
              for (k in seq_len(length(tmpind))) {
                lines(tmp[tmpind[k],], col = colo[jj])
              }
            }
          }
        }
      }
    }
    if (!is.na(filename))
      dev.off()
  }

#' Plotting results from estimating the cluster number
#'
#' This function visualizes the output from estimClustNumber, and there 
#' particularly the
#' two validity indices Minimum Centroid Distance and Xie Beni Index.
#'
#' @param ClustInd Matrix with values from validity indices
#' @return Multiple panels showing expression profiles of clustered features 
#' passing the minMem threshold
#' @examples
#' data("artificial_clusters")
#' dat <- averageCond(artificial_clusters, 5, 10)
#' dat <- scale(dat)
#' dat <- cbind(dat, 1)
#' ClustInd <- estimClustNum(dat, 6)
#' estimClust.plot(ClustInd)
#' @export
#' @references
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering 
#' of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 
#' 10.1093/bioinformatics/bty224. PMID: 29635359.
#'
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and 
#' the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 
#' '2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#'
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters 
#' for fuzzy c-means cluster analysis. Bioinformatics. 2010 
#' Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. 
#' PMID: 20880957.
estimClust.plot <- function(ClustInd) {
  par(mfrow = c(1, 3))
  maxClust <- nrow(ClustInd) + 2
  plot(
    seq(3,maxClust,1),
    ClustInd[seq_len(nrow(ClustInd)), "MinCentroidDist_VSClust"],
    col = 2 ,
    type = "b",
    main = "Min. centroid distance\n(Highest jump is best)",
    xlab = "Number of clusters",
    ylab = "Index",
    ylim = c(min(ClustInd[, grep("MinCentroidDist",
                                 colnames(ClustInd))], na.rm = TRUE),
             max(ClustInd[, grep("MinCentroidDist",
                                 colnames(ClustInd))], na.rm = TRUE))
  )
  lines(seq(3,maxClust,1), ClustInd[seq_len(nrow(ClustInd)), "MinCentroidDist_FCM"],
        col = 3, type = "b")
  dmindist <- optimalClustNum(ClustInd)
  points(dmindist,
         ClustInd[dmindist - 2, "MinCentroidDist_VSClust"],
         pch = 15,
         col = 1,
         cex = 2)
  legend(
    "topright",
    legend = c("VSClust", "Standard"),
    lty = c(1, 1),
    col = 2:3
  )
  grid(NULL, NA, lwd = 1, col = 1)
  plot(
    seq(3,maxClust,1),
    ClustInd[seq_len(nrow(ClustInd)), "XieBeni_VSClust"],
    col = 2,
    type = "b",
    main = "Xie-Beni index\n(Lowest is best)",
    xlab = "Number of clusters",
    ylab = "Index",
    ylim = c(min(ClustInd[, grep("XieBeni", colnames(ClustInd))], na.rm =
                   TRUE),
             max(ClustInd[, grep("XieBeni", colnames(ClustInd))], na.rm =
                   TRUE))
  )
  lines(seq(3,maxClust,1), ClustInd[seq_len(nrow(ClustInd)), "XieBeni_FCM"], type =
          "b", col = 3)
  dxiebeni <- optimalClustNum(ClustInd, index = "XieBeni")
  points(dxiebeni,
         ClustInd[dxiebeni - 2, "XieBeni_VSClust"],
         pch = 15,
         col = 1,
         cex = 2)
  legend(
    "topright",
    legend = c("VSClust", "Standard"),
    lty = c(1, 1),
    col = 2:3
  )
  grid(NULL, NA, lwd = 1, col = 1)
  plot(
    3:maxClust,
    ClustInd[seq_len(nrow(ClustInd)), "NumVSClust"],
    col = 2,
    type = "b",
    main = "Total number of assigned features",
    xlab = "Number of clusters",
    ylab = "Assigned features",
    ylim = c(min(ClustInd[, grep("Num", colnames(ClustInd))], na.rm =
                   TRUE),
             max(ClustInd[, grep("Num", colnames(ClustInd))], na.rm = TRUE))
  )
  lines(3:maxClust, ClustInd[seq_len(nrow(ClustInd)), "NumFCM"], type =
          "b", col = 3)
  legend(
    "topright",
    legend = c("VSClust", "Standard"),
    lty = c(1, 1),
    col = 2:3
  )
  # finally plot
  p <- recordPlot()
}
