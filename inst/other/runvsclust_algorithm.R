library(Rcpp)
# basic algorithm
sourceCpp("source/c++means_missing_values.cpp")
# function to run c++ code
source("vsclust_algorithm.R")

# test function to run vsclust with one fuzzifier and multiple ones
# parameters are the number of features and the fracation of missing values
# The data set contains 2 well-separated clusters
run_test <- function(n=1000, na=0, ...) {
  
  dat <- matrix(rnorm(n)+rep(c(-1,1),n/2), n/10)
  dat[sample(1:length(dat), size=na*length(dat))] <- NA
  par(mfrow=c(1,2))
  # setting arbitrary fuzzifier values
  m <- 1.1 + runif(n/10)
  
  cmeans_out <- vsclust_algorithm(dat,centers=2, m=1.5, ...)
  hist(cmeans_out$membership, 100, xlab="All membership values", col="#333333", border=NA, 
       main="Fuzzifier set to 1.5")
  cmeans_out2 <- vsclust_algorithm(dat,centers=2, m=m, ...)
  hist(cmeans_out2$membership, 100, xlab="All membership values", col="#333333", border=NA, 
       main="Individual fuzzifiers")
  par(mfrow=c(1,1))
  return(list(cmeans1=cmeans_out, cmeans2=cmeans_out2))
}

# run test
cmp_data <- run_test()



###### other stuff ###

# Checking run time by increasing feature number
t <- NULL
for (i in rep(10^(seq(2,7,1)), each=5)) {
  t <- append(t,system.time(run_test(i, na=0.0))[3])
}
plot(t)

# If available, running the old code from the modified e1071 library
run_test2 <- function(n=1000, ...) {
  
  dat <- matrix(rnorm(n)+ rep(c(-1,1),n/2), n/10)
  #dat <- t(scale(t(dat)))
  
  #pca_out <- princomp(dat)
  #plot(pca_out$scores)
  par(mfrow=c(1,2))
  m <- 11 + runif(n/10)
  
  cmeans_out <- cmeans(dat,centers=2, m=1.5, ...)
  hist(as.vector(cmeans_out$membership))
  
  cmeans_out <- cmeans(dat,centers=2, m=m, ...)
  hist(as.vector(cmeans_out$membership))
}

t <- NULL
for (i in rep(10^(seq(2,7,1)), each=5)) {
  t <- append(t,system.time(run_test2(i))[3])
}
plot(t)
