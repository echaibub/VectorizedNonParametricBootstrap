##############################################################
## Additional functions implementing a vectorized version
## of the approximate bootstrap based on Poisson frequencies
##############################################################

## Computes bootstrap weights from Poison counts
## N: number of data points
## B: number of bootstrap replications
## returns B by N matrix of bootstrap weights
##
BootPoissonWeights <- function(N, B) {
  counts <- rpois(N * B, lambda = 1)
  counts <- matrix(counts, B, N)
  aux <- apply(counts, 1, sum)
  counts/aux
}


## Same as SampleBootCor, except that W is B by N,
## instead of N by B.
##
## Computes the correlation coefficient in vectorized form
## W: bootstrap weights matrix
## x1: variable 1 data vector
## x2: variable 2 data vector
## returns a vector of length B of correlation coefficients
##
SampleBootCor2 <- function(W, x1, x2) {
  xbar.1 <- W %*% x1
  xbar.2 <- W %*% x2
  s2.1 <- x1^2
  s2.1 <- W %*% s2.1
  s2.1 <- s2.1 - xbar.1^2
  s2.2 <- x2^2
  s2.2 <- W %*% s2.2
  s2.2 <- s2.2 - xbar.2^2  
  s.12 <- x1 * x2
  s.12 <- W %*% s.12
  s.12 <- s.12 - xbar.1 * xbar.2
  out <- s.12/sqrt(s2.1 * s2.2)
}


## Runs the vectorized approximate bootstrap based on Poisson frequencies 
## N: number of data points
## B: number of bootstrap replications
## theta: function of the bootstrap weights, computing the statistic 
##        of interest (in vectorized form)
## ... : additional arguments passed to the theta function 
## returns a vector of bootstrap replications
##
VectorizedPoissonBootstrap <- function(N, B, theta, ...) {
  call <- match.call()
  W <- BootPoissonWeights(N, B)
  theta.star <- theta(W, ...)
  as.vector(theta.star)
}

