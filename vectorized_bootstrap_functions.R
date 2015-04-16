#######################################
## basic R functions
#######################################

library(bootstrap)
library(boot)

## Computes bootstrap weights from multinomial counts
## N: number of data points
## B: number of bootstrap replications
## returns N by B matrix of bootstrap weights
##
BootWeights <- function(N, B) {
  counts <- rmultinom(B, N, rep(1/N, N))
  counts/N
}

## Runs the vectorized bootstrap
## N: number of data points
## B: number of bootstrap replications
## theta: function of the bootstrap weights, computing the statistic 
##        of interest (in vectorized form)
## ... : additional arguments passed to the theta function 
## returns a vector of bootstrap replications
##
VectorizedBootstrap <- function(N, B, theta, ...) {
  call <- match.call()
  W <- BootWeights(N, B)
  theta.star <- theta(W, ...)
  as.vector(theta.star)
}

## Computes the correlation coefficient in vectorized form
## W: bootstrap weights matrix
## x1: variable 1 data vector
## x2: variable 2 data vector
## returns a vector of length B of correlation coefficients
##
SampleBootCor <- function(W, x1, x2) {
  xbar.1 <- crossprod(x1, W)
  xbar.2 <- crossprod(x2, W)
  s2.1 <- x1^2
  s2.1 <- crossprod(s2.1, W)
  s2.1 <- s2.1 - xbar.1^2
  s2.2 <- x2^2
  s2.2 <- crossprod(s2.2, W)
  s2.2 <- s2.2 - xbar.2^2  
  s.12 <- x1 * x2
  s.12 <- crossprod(s.12, W)
  s.12 <- s.12 - xbar.1 * xbar.2
  s.12/sqrt(s2.1 * s2.2) 
}

## Correlation coefficient function for the bootstrap R function call
## idx: re-sampled data index
## x1: variable 1 data vector
## x2: variable 2 data vector
## returns the correlation between the re-sampled versions of x1 and x2
##
ThetaCor <- function(idx, x1, x2) {
  cor(x1[idx], x2[idx])
}

## Correlation coefficient function for the boot R function call
## x: data matrix containing variables 1 and 2 
## idx: re-sampled data index
## returns the correlation between the re-sampled versions of x1 and x2
##
ThetaCor2 <- function(x, idx) {
  cor(x[idx, 1], x[idx, 2])
}

## Loop implementation for bootstrapping Pearson's correlation
## B: number of bootstraps
## x1: variable 1 data vector
## x2: variable 2 data vector
## returns a vector of length B of correlation coefficients
##
LoopCorBootstrap <- function(B, x1, x2) {
  N <- length(x1)
  out <- rep(NA, B) 
  for (i in seq(B)) {
    idx <- sample(N, replace = TRUE)
    out[i] <- cor(x1[idx], x2[idx])
  }
  out
}

## Similar to VectorizedBootstrap(), but instead of generating W internaly
## it takes W as an argument. (This function was used for the analysis in
## Fig. 4)
##
## Runs the vectorized bootstrap
## W: bootstrap weight matrix
## theta: function of the bootstrap weights, computing the statistic 
##        of interest (in vectorized form)
## ... : additional arguments passed to the theta function 
## returns a vector of bootstrap replications
##
VectorizedBootstrap2 <- function(W, theta, ...) {
  call <- match.call()
  theta.star <- theta(W, ...)
  as.vector(theta.star)
}
