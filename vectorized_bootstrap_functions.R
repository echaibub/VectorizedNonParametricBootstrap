###############################################################
## basic R functions used for the comparisons in the main text
###############################################################

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

## Runs the vectorized bootstrap on two sample problems
## N1: number of data points on group 1
## N2: number of data points on group 2
## B: number of bootstrap replications
## theta: function of the bootstrap weights, computing the statistic 
##        of interest (in vectorized form)
## ... : additional arguments passed to the theta function 
## returns a vector of bootstrap replications
##
Vectorized2SampleBootstrap <- function(N1, N2, B, theta, ...) {
  call <- match.call()
  W1 <- BootWeights(N1, B)
  W2 <- BootWeights(N2, B)
  theta.star <- theta(W1, W2, ...)
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


## Generates the results for Fig. 4
##
ReplicateVecTiming <- function(K = 3, N.seq, x1, x2) {
  n.N <- length(N.seq)
  at1.W <- at2.W <- at3.W <- matrix(0, n.N, 3)
  at1.M <- at2.M <- at3.M <- matrix(0, n.N, 3)  
  for (k in seq(K)) {
    t1.W <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
    t1.M <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
    t2.W <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
    t2.M <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
    t3.W <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
    t3.M <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
    for(i in seq(n.N)) {
      cat("k, i = ", c(k, i), "\n")
      x1 <- X1[seq(N.seq[i])]
      x2 <- X2[seq(N.seq[i])] 
      cat("1e+4, ", "\n")
      W1 <- system.time(W <- BootWeights(N.seq[i], 10000))
      M1 <- system.time(out <- VectorizedBootstrap2(W, SampleBootCor, x1, x2))
      t1.W[i,] <- W1[1:3]
      t1.M[i,] <- M1[1:3]
      cat("1e+5, ", "\n")
      W2 <- system.time(W <- BootWeights(N.seq[i], 100000))
      M2 <- system.time(out <- VectorizedBootstrap2(W, SampleBootCor, x1, x2))
      t2.W[i,] <- W2[1:3]
      t2.M[i,] <- M2[1:3]
      cat("1e+6, ", "\n")
      W3 <- system.time(W <- BootWeights(N.seq[i], 1000000))
      M3 <- system.time(out <- VectorizedBootstrap2(W, SampleBootCor, x1, x2))
      t3.W[i,] <- W3[1:3]
      t3.M[i,] <- M3[1:3]
    }
    at1.W <- at1.W + t1.W 
    at2.W <- at2.W + t2.W
    at3.W <- at3.W + t3.W
    at1.M <- at1.M + t1.M 
    at2.M <- at2.M + t2.M
    at3.M <- at3.M + t3.M
  }
  list(at1.W = at1.W/K, at2.W = at2.W/K, at3.W = at3.W/K,
       at1.M = at1.M/K, at2.M = at2.M/K, at3.M = at3.M/K)
}


##############################################
## Additional function used in the tutorial
##############################################

## mean example

ThetaMean <- function(idx, x) {
  mean(x[idx])
}

ThetaMean2 <- function(x, idx) {
  mean(x[idx])
}

LoopMeanBootstrap <- function(B, x) {
  N <- length(x)
  out <- rep(NA, B) 
  for (i in seq(B)) {
    idx <- sample(N, replace = TRUE)
    out[i] <- mean(x[idx])
  }
  out
}

SampleBootMean <- function(W, x) {
  crossprod(x, W)
}


## Welch's statistic example

LoopWelchsBootstrap <- function(B, x1, x2) {
  n1 <- length(x1)
  n2 <- length(x2)
  n <- n1 + n2
  stat <- rep(NA, B) 
  x.bar <- mean(c(x1, x2))
  x1.tilda <- x1 - mean(x1) + x.bar
  x2.tilda <- x2 - mean(x2) + x.bar
  for (i in seq(B)) {
    i1 <- sample(n1, replace = TRUE)
    i2 <- sample(n2, replace = TRUE)
    stat[i] <- WelchsTestStat(x1.tilda[i1], x2.tilda[i2])
  }
  stat
}

WelchsTestStat <- function(x1, x2) {
  n1 <- length(x1)
  n2 <- length(x2)
  s2.1 <- var(x1)
  s2.2 <- var(x2)
  xbar.1 <- mean(x1)
  xbar.2 <- mean(x2)
  s.12 <- sqrt((s2.1/n1) + (s2.2/n2))
  (xbar.1 - xbar.2)/s.12
}

ThetaWelchs <- function(x, f) {
  i1 <- which(x[, 2] == 1)
  i2 <- which(x[, 2] == 2)
  n1 <- length(i1)
  n2 <- length(i2)
  x.tilda <- x
  x.tilda[i1, 1] <- x[i1, 1] - mean(x[i1, 1]) + mean(x[, 1])
  x.tilda[i2, 1] <- x[i2, 1] - mean(x[i2, 1]) + mean(x[, 1])
  x <- x.tilda
  xbar.1 <- sum(x[i1, 1] * f[i1])/sum(f[i1])
  s2.1 <- sum(x[i1, 1]^2 * f[i1])/sum(f[i1]) - xbar.1^2
  s2.1 <- n1 * s2.1/(n1 - 1)
  xbar.2 <- sum(x[i2, 1] * f[i2])/sum(f[i2])
  s2.2 <- sum(x[i2, 1]^2 * f[i2])/sum(f[i2]) - xbar.2^2  
  s2.2 <- n2 * s2.2/(n2 - 1)
  s.12 <- sqrt((s2.1/n1) + (s2.2/n2))
  (xbar.1 - xbar.2)/s.12
}

SampleBootWelchs <- function(W1, W2, x1, x2) {
  n1 <- length(x1)
  n2 <- length(x2)
  x.bar <- mean(c(x1, x2))
  x1.tilda <- x1 - mean(x1) + x.bar
  x2.tilda <- x2 - mean(x2) + x.bar
  xbar.1 <- crossprod(x1.tilda, W1)
  xbar.2 <- crossprod(x2.tilda, W2)  
  s2.1 <- x1.tilda^2
  s2.1 <- crossprod(s2.1, W1)
  s2.1 <- s2.1 - xbar.1^2 
  s2.1 <- n1 * s2.1/(n1 - 1)
  s2.2 <- x2.tilda^2
  s2.2 <- crossprod(s2.2, W2)
  s2.2 <- s2.2 - xbar.2^2 
  s2.2 <- n2 * s2.2/(n2 - 1)
  s.12 <- sqrt((s2.1/n1) + (s2.2/n2))
  (xbar.1 - xbar.2)/s.12
}





