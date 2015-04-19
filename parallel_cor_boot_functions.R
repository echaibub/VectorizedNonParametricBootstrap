########################################################
## parallel functions for bootstrapping Pearson's
## correlation coefficient, leveraging the parallel R 
## package standard capabilities
## #####################################################


library("parallel")


## Implements the parLapply parallel version of the 
## "for loop" computation of Pearson's correlation coefficient.
## B: number of bootstrap replications
## x1: variable 1 data vector
## x2: variable 2 data vector
## mc: number of cores
## returns a vector of length B of correlation coefficients
##
ParallelLoopCorBootstrap <- function(B, x1, x2, mc = 4) {
  cl <- makeCluster(mc, type = "PSOCK")
  results <- parLapply(cl, SplitReplications(B, mc), fun = LoopCorBootstrap, x1, x2) 
  results <- unlist(results)
  stopCluster(cl)
  results
}


## Implements the mclapply parallel version of the 
## "for loop" computation of Pearson's correlation coefficient.
## B: number of bootstrap replications
## x1: variable 1 data vector
## x2: variable 2 data vector
## mc: number of cores
## returns a vector of length B of correlation coefficients
##
ParallelLoopCorBootstrap2 <- function(B, x1, x2, mc = 4) {
  results <- mclapply(SplitReplications(B, mc), LoopCorBootstrap, x1, x2) 
  unlist(results)
}


## Implements the parLapply parallel version of the 
## vectorized computation of Pearson's correlation coefficient.
## N: number of data points
## B: number of bootstrap replications
## x1: variable 1 data vector
## x2: variable 2 data vector
## mc: number of cores
## returns a vector of length B of correlation coefficients
##
ParallelVectorizedCorBootstrap <- function(B, N, x1, x2, mc = 4) {
  cl <- makeCluster(mc, type = "PSOCK")
  results <- parLapply(cl, SplitReplications(B, mc), fun = VectorizedCorBootstrap, N, x1, x2) 
  results <- unlist(results)
  stopCluster(cl)
  results
}


## Implements the mclapply parallel version of the 
## vectorized computation of Pearson's correlation coefficient.
## N: number of data points
## B: number of bootstrap replications
## x1: variable 1 data vector
## x2: variable 2 data vector
## mc: number of cores
## returns a vector of length B of correlation coefficients
##
ParallelVectorizedCorBootstrap2 <- function(B, N, x1, x2, mc = 4) {
  results <- mclapply(SplitReplications(B, mc), VectorizedCorBootstrap, N, x1, x2) 
  unlist(results)
}


## Splits the number of bootstrap replications, B, into k roughly 
## equal sized groups.
## returns a vector of length k with the group sizes
##
SplitReplications <- function(x, k) {
  aux1 <- x %/% k
  aux2 <- x %% k
  c(rep(aux1, k - 1), aux1 + aux2)
}

