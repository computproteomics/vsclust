#' Synthetic/artificial data comprising 5 clusters
#'
#' 10-dimensional data set with 500 simulating features measured over 5 replicates each, comprising a total of 50 samples. 
#' The first 250 features were modeled through normal distributions shifted in the 10-dimensional space to form 5 different clusters.
#' The 2nd half of the features were modeled through a normal distribution around the origin and thus should be assigned to any cluster
#'
#' @docType data
#' @usage data(artificial_clusters)
#' @format A data frame consisting of 500 features distributed over 5 clusters and being replicated 5 times each
#' @source Protein Research Group, University of Southern Denmark, Odense
"artificial_clusters"
