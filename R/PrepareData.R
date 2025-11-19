#' Paired statistical testing
#'
#' Statistical testing and variance estimation in multi-dimensional data set. 
#' given by a matrix. This functions runs LIMMA paired tests and
#' calculated the shrunken variance estimates.
#'
#' @param Data a numeric data matrix with columns as samples. Different 
#' experimental conditions are grouped together in their replicates. The number 
#' of samples per group needs to be identical
#' @param NumCond Number of different experimental conditions
#' @param NumReps Number of replicates per experimental condition
#' @return List containing the objects
#' @return `qvalues` false discovery rates after correction for multiple testing 
#' (`qvalue` method from `qvalue` library)
#' @return `Sds` General standard deviation within replicates after using 
#' shrinkage (eBayes) by LIMMA
#' @examples
#' #' # Generate some random data with three different experimental conditions
#' data <- matrix(rnorm(seq_len(1500)), nrow=100)
#' # Run statistical testing
#' stat_out <- SignAnalysisPaired(data, 3, 5)
#' # Histogram of qvalues comparing the second to the first condition
#' hist(stat_out$qvalues[,1], 50, xlab="q-values")
#' @import limma
#' @import stats
#' @importFrom qvalue qvalue
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
#' @export
SignAnalysisPaired <- function(Data, NumCond, NumReps) {
    ##########################################################
    # significance analysis
    MAData <- Data[, 2:(NumCond)] - Data[, 1]
    for (i in seq_len(NumReps - 1))
        MAData <-
            cbind(MAData, Data[, (i * NumCond + 1) + seq_len(NumCond - 1)] - 
                      Data[, (i * NumCond + 1)])
    rownames(MAData) <- rownames(Data)
    MAReps <- rep(seq_len(NumCond - 1), NumReps)
    if (is.null(rownames(MAData)))
        rownames(MAData) <- paste0("feature", seq_len(nrow(MAData)))
    ##limma with ratios
    design <- plvalues <- NULL
    for (c in (seq_len(NumCond - 1))) {
        design <- cbind(design, as.numeric(MAReps == c))
    }
    lm.fittedMA <- lmFit(MAData, design)
    lm.bayesMA <- eBayes(lm.fittedMA)
    topTable(lm.bayesMA)
    plvalues <- lm.bayesMA$p.value
    qvalues <-
        matrix(
            NA,
            nrow = nrow(plvalues),
            ncol = ncol(plvalues),
            dimnames = dimnames(plvalues)
        )
    # qvalue correction
    for (i in seq_len(ncol(plvalues))) {
        tqs <- tryCatch(
            qvalue(na.omit(plvalues[, i]))$qvalues,
            error = function(e)
                NULL
        )
        if (length(tqs) > 0) {
            qvalues[names(tqs), i] <- tqs
        }
        else {
            qvalues[names(tqs), i] <- NA
        }
    }
    
    return(list(qvalues = qvalues, Sds = sqrt(lm.bayesMA$s2.post)))
}

#' Unpaired statistical testing
#'
#' Statistical testing and variance estimation in multi-dimensional data set. 
#' given by a matrix. This functions runs LIMMA paired tests and
#' calculated the shrunken variance estimates.
#'
#' @param Data a numeric data matrix with columns as samples. Different 
#' experimental conditions are grouped together in their replicates. The number 
#' of samples per group needs to be identical
#' @param NumCond Number of different experimental conditions
#' @param NumReps Number of replicates per experimental condition
#' @return List containing the objects
#' @return `pvalues` p-values before correction for multiple testing
#' @return `qvalues` false discovery rates after correction for multiple testing 
#' (`qvalue` method from `qvalue` library)
#' @return `Sds` General standard deviation within replicates after using 
#' shrinkage by LIMMA
#' @examples
#' #' # Generate some random data
#' data <- matrix(rnorm(seq_len(1000)), nrow=100)
#' # Run statistical testing
#' stat_out <- SignAnalysis(data, 2, 5)
#' # Histogram of qvalues (no significant events)
#' hist(stat_out$qvalues, 50, xlab="q-values")
#' @import limma
#' @import stats
#' @importFrom qvalue qvalue
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
SignAnalysis <- function(Data, NumCond, NumReps) {
    ##########################################################
    if (is.null(rownames(Data)))
        rownames(Data) <- paste0("feature", seq_len(nrow(Data)))
    # significance analysis
    Reps <- rep(seq_len(NumCond), NumReps)
    design <- model.matrix( ~ 0 + factor(Reps - 1))
    colnames(design) <- paste("i", c(seq_len(NumCond)), sep = "")
    contrasts <- NULL
    First <- 1
    for (i in (seq_len(NumCond))[-First])
        contrasts <- append(contrasts,
                            paste(colnames(design)[i], "-", colnames(design)[First], 
                                  sep = ""))
    contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)
    lm.fitted <- lmFit(Data, design)
    
    lm.contr <- contrasts.fit(lm.fitted, contrast.matrix)
    lm.bayes <- eBayes(lm.contr)
    #topTable(lm.bayes)
    plvalues <- lm.bayes$p.value
    qvalues <- matrix(
        NA,
        nrow = nrow(plvalues),
        ncol = ncol(plvalues),
        dimnames = dimnames(plvalues)
    )
    # qvalue correction
    for (i in seq_len(ncol(plvalues))) {
        tqs <- tryCatch(
            qvalue(na.omit(plvalues[, i]))$qvalues,
            error = function(e)
                NULL
        )
        # print(tqs)
        if (length(tqs) > 0) {
            qvalues[names(tqs), i] <- tqs
        }
        else {
            qvalues[names(tqs), i] <- NA
        }
    }
    return(list(
        pvalues = plvalues,
        qvalues = qvalues,
        Sds = sqrt(lm.bayes$s2.post)
    ))
}

#' Visualize using principal component analysis (both loadings and scoring) 
#' including the variance from the replicates
#'
#' The loading plot shows all features and their scaled variance. This provides 
#' an idea of the intrinsic noise in the data.
#'
#' @param data Matrix of data frame with numerical values. Columns corresponds 
#' to samples
#' @param NumReps Number of replicates per experimental condition
#' @param NumCond Number of different experimental conditions
#' @param Sds Standard deviation for each features. Usually using the one from 
#' LIMMA
#' @return Loading and scoring plots that include feature variance
#' @examples
#' data <- matrix(rnorm(1000), nrow=100)
#' pcaWithVar(data, NumCond=2, NumReps=5, Sds=1)
#' @import stats
#' @import graphics
#' @importFrom shiny need isRunning validate
#' @importFrom grDevices rainbow
#' @export
#' @references
#' Schwaemmle V, Jensen ON. VSClust: feature-based variance-sensitive clustering 
#' of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. 
#' doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
#'
#' Schwaemmle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and 
#' the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 
#' 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.
#'
#' Schwaemmle V, Jensen ON. A simple and fast method to determine the parameters 
#' for fuzzy c-means cluster analysis. Bioinformatics. 2010 
#' Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29.
#' PMID: 20880957.
pcaWithVar <- function(data, NumReps, NumCond, Sds = 1) {
    # Remove columns with only 5% of the data
    plab <- rep(seq_len(NumCond), NumReps)
    plab <- plab[colSums(!is.na(data)) > 0.05 * nrow(data)]
    pcaDat <- data
    if (ncol(pcaDat) != NumCond * NumReps)
        stop("Wrong number of conditions and/or replicates!")
    pcaDat <- data[, colSums(!is.na(data)) > 0.05 * nrow(data)]
    pcaDat <- (pcaDat[complete.cases(pcaDat), ])
    if (shiny::isRunning()){
        validate(need(
            length(pcaDat) > 0,
            "Principal component analysis not shown as too many missing values"
        ))
        validate(need(
            nrow(pcaDat) > 10,
            "Principal component analysis not shown as too many missing values"
        ))
    } else {
        if (length(pcaDat) == 0 || nrow(pcaDat) < 10){
            warning("Principal component analysis not possible due to too many missing values")
            return(NULL)
        }
    }
    pca <- prcomp(pcaDat, scale = TRUE, retx = TRUE)
    # scores <- pca$x
    # loadings <- pca$rotation
    scores <- pca$rotation
    loadings <- pca$x
    par(mfrow = c(1, 2))
    # Set missing values to maximum std. dev.
    Sds[is.na(Sds)] <- max(Sds, na.rm = TRUE)
    # Scaling suitable for visualization
    Sds <- sqrt(Sds)
    plot(
        loadings,
        cex = Sds,
        pch = 16,
        col = paste("#000000", sprintf("%02X", as.integer(
            255 / max(1 / Sds) / Sds
        )), sep = "")
    )
    title(main = "Principal component analysis of data set (loadings)", sub =
              "The point size corresponds to the estimated standard deviation")
    plot(scores, pch = 19, col = rainbow(NumCond)[plab])
    title(main = "Principal component analysis of data set (scores)", sub =
              "Colors denote different conditions")
    legend(
        "topright",
        paste("Condition", seq_len(NumCond)),
        pch = rep(19, NumCond),
        col = rainbow(NumCond)[seq_len(NumCond)]
    )
}

#' Determine optimal cluster number from validity index
#'
#' Calculated the optimal number from expected behavior of the indices. This 
#' would be a large decay for the Minimum Centroid Distance and a minimum for 
#' the Xie Beni index
#'
#' @param ClustInd Output from estimClustNum providing the calculated cluster 
#' validity indices
#' @param index Either "MinCentroidDist" or "XieBeni"
#' @param method Either "VSClust" or "FCM" for standard fuzzy c-means clustering
#' @return optimal cluster number
#' @examples
#'   data("artificial_clusters")
#'   dat <- averageCond(artificial_clusters, 5, 10)
#'   dat <- scale(dat)
#' dat <- cbind(dat, 1)
#' ClustInd <- estimClustNum(dat, 6)
#' optimalClustNum
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
#' for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. 
#' doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
optimalClustNum <-
    function(ClustInd,
             index = "MinCentroidDist",
             method = "VSClust") {
        allowedInd <- c("XieBeni", "MinCentroidDist")
        allowedMethod <- c("FCM", "VSClust")
        if (!any(index == allowedInd)) {
            stop(paste("index needs to be one of", paste(allowedInd, collapse = " ")))
        }
        if (!any(method == allowedMethod)) {
            stop(paste(
                "method needs to be one of",
                paste(allowedMethod, collapse = " ")
            ))
        }
        
        tClustInd <- ClustInd[, grep(index, colnames(ClustInd))]
        tClustInd <- tClustInd[, grep(method, colnames(tClustInd))]
        opt_val <- NULL
        if (length(tClustInd) < 3)
            stop("Minimal length of ClustInd vector is 3")
        if (index == "MinCentroidDist") {
            opt_val <-
                which.max(tClustInd[seq_len(length(tClustInd) - 1)] - 
                              tClustInd[2:length(tClustInd)])
        } else if (index == "XieBeni") {
            opt_val <- which.min(tClustInd[seq_len(length(tClustInd))])
            
        }
        return(as.numeric(opt_val + 2))
    }

