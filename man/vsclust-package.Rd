\name{vsclust-package}
\alias{vsclust-package}
\alias{vsclust}
\docType{package}
\title{
  VSClust provides a powerful method to run variance-sensitive clustering
}
\description{
  Clustering of high-dimensional quantitative data with data points that come with multiple measurements.
  In this clustering method, each feature is represented by a) its 
  quantitative profile and b) its variance. Hence, the incertainty about a measurement enter
  in the determination of the most common patterns. This methods is both insensitive to 
  noisy measurements and avoids finding clusters in homogeneously distributed data.
}
\details{
  The functions in this package comprise (i) methods to prepare the data for cluster analysis like statistical analysis (`SignAnal`and SignPairedAnal`), PCA (`PCAwithVar`), (ii) direct application of the clustering algorithm on a (standardized) data matrix (`vsclust_algorithm`), (iii) for the further evaluation and visualization (such as `calcBHI` and `mfuzz.plot`), and (iv) wrappers for the over workflow including statistical preparation (`statWrapper`), estimation of the cluster number (`estimClustNum`), running the clustering (`runClustWrapper`) and functional evaluation (`runFuncEnrich`).
}
\author{
Maintainer: Veit Schw\"ammle" <veits@bmb.sdu.dk>
}
\references{
  - Schw\"ammle V, Hagensen CE. A Tutorial for Variance-Sensitive Clustering and the Quantitative Analysis of Protein Complexes. Methods Mol Biol. 2021;2228:433-451. doi: 10.1007/978-1-0716-1024-4_30. PMID: 33950508.

  - Schw\"ammle V, Jensen ON. VSClust: feature-based variance-sensitive clustering of omics data. Bioinformatics. 2018 Sep 1;34(17):2965-2972. doi: 10.1093/bioinformatics/bty224. PMID: 29635359.
  
  - Schw\"ammle V, Jensen ON. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. Bioinformatics. 2010 Nov 15;26(22):2841-8. doi: 10.1093/bioinformatics/btq534. Epub 2010 Sep 29. PMID: 20880957.
}
\keyword{ vsclust }
\seealso{

}
\examples{
}
