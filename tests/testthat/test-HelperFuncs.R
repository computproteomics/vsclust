# Both tests check whether most of the clustered feature have been found (50% of the 500 features)
test_that("vsclust_algorithm", {
  data("artificial_clusters")
  dat <- averageCond(artificial_clusters, 5, 10)
  dat <- scale(dat)

  clust_out <- vsclust_algorithm(dat, centers = 6, m = 1.55)
  expect_equal(    sum(apply(clust_out$membership, 1, max) > 0.5), 200, tolerance = 10)
})


test_that("clust_comp", {
  data("artificial_clusters")
  dat <- averageCond(artificial_clusters, 5, 10)
  clust_out <- ClustComp(dat, NClust = 6, Sds = 1)
  sum(apply(clust_out$Bestcl$membership, 1, max) > 0.5)
  expect_equal(  sum(apply(clust_out$Bestcl$membership, 1, max) > 0.5), 214, tolerance = 10)
})

test_that("sign_analysis_paired", {
  set.seed(0)
  data <- matrix(rnorm(1:1500), nrow=100)
  # make artificial regulations
  data[1:20,c(2,5,8,11,14)] <- data[1:20, c(1,4,7,10,13)] + 3
  data[21:40,c(2,5,8,11,14)] <- data[21:40, c(1,4,7,10,13)] - 3
  # Run statistical testing
  stat_out <- SignAnalysisPaired(data, 3, 5)
  # Histogram of qvalues comparing the second to the first condition
  expect_equal(sum(stat_out$qvalues[,1] < 0.01), 42)
})

test_that("sign_analysis_unpaired", {
  set.seed(1)
  data <- matrix(rnorm(1:1500), nrow=100)
  # make artificial regulations
  data[1:20,c(2,5,8,11,14)] <- data[1:20,c(2,5,8,11,14)] + 3
  data[21:40,c(2,5,8,11,14)] <- data[21:40,c(2,5,8,11,14)] - 3
  # Run statistical testing
  stat_out <- SignAnalysis(data, 3, 5)
  # Histogram of qvalues comparing the second to the first condition
  expect_equal(sum(stat_out$qvalues[,1] < 0.01), 40)
})

test_that("prepare_for_vsclust", {
  set.seed(1)
  data <- matrix(rnorm(1:1500), nrow=100)
  # make artificial regulations
  data[1:20,c(2,5,8,11,14)] <- data[1:20,c(2,5,8,11,14)] + 3
  data[21:40,c(2,5,8,11,14)] <- data[21:40,c(2,5,8,11,14)] - 3
  # Run statistical testing
  prepared <- PrepareForVSClust(data, 5,3, isStat=TRUE)
  expect_equal(sum(prepared$statFileOut[,"qvalue BvsA"] < 0.01), 40)
})

test_that("estimate_clust_num", {
  set.seed(0)
  data("artificial_clusters")
  dat <- averageCond(artificial_clusters, 5, 10)
  dat <- scale(dat)
  dat <- cbind(dat, 1)
  ClustInd <- estimClustNum(dat, 10)
  expect_equal(as.numeric(optimalClustNum(ClustInd, index="MinCentroidDist", method="VSClust")), 6)
})

