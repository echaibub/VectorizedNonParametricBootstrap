
#######################################################
## Script generating the results in Section 3 of the 
## manuscript
#######################################################

library(devtools)
source_url("https://raw.githubusercontent.com/echaibub/VectorizedNonParametricBootstrap/master/vectorized_bootstrap_functions.R")

#######################################
## American law schools example
#######################################

## load the American law schools data from the bootstrap R package
##
data(law82) 
N <- nrow(law82)
x1 <- law82$LSAT
x2 <- law82$GPA
x <- cbind(x1, x2)

## just check the data-resampling and vectorized implementations 
## on 100,000 bootstraps
##
B <- 100000
system.time(standard1 <- bootstrap(seq(N), B, ThetaCor, x1, x2)$thetastar)
system.time(standard2 <- boot(x, ThetaCor2, B)$t[, 1])
system.time(vectorized <- VectorizedBootstrap(N, B, SampleBootCor, x1, x2))
par(mfrow = c(1, 3))
hist(standard1, nclass = 50, xlim = c(0.4, 1))
hist(standard2, nclass = 50, xlim = c(0.4, 1))
hist(vectorized, nclass = 50, xlim = c(0.4, 1))
par(mfrow = c(1, 1))

## Run the comparison of the vectorized implementation against the 
## data re-sampling implementations in: (i) the bootstrap R package 
## (ii) boot R package; and (iii) on a strait foward for loop 
## implementation.
##
## The implementations are compared over a grid of 30 distinct number 
## of bootstraps values ranging from 1,000 to 1,000,000.
##
B.seq <- round(exp(seq(log(1e+3), log(1e+6), length.out = 30)), 0)
set.seed(123456789)
n.B <- length(B.seq)
t0.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t1.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t2.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t3.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
for(i in seq(n.B)) {
  cat("i = ", i, "\n")
  t0 <- system.time(out0 <- bootstrap(seq(N), B.seq[i], ThetaCor, x1, x2)$thetastar)
  t1 <- system.time(out1 <- VectorizedBootstrap(N, B.seq[i], SampleBootCor, x1, x2))
  t2 <- system.time(out2 <- LoopCorBootstrap(B.seq[i], x1, x2))
  t3 <- system.time(out3 <- boot(x, ThetaCor2, B.seq[i])$t[, 1])
  t0.time[i,] <- t0[1:3]
  t1.time[i,] <- t1[1:3]
  t2.time[i,] <- t2[1:3]
  t3.time[i,] <- t3[1:3]
}

## Generate Figure 1 of the manuscript.
##
cl <- 0.9
ca <- 0.9
cm <- 0.85
cp <- 0.75
k <- 1
nc <- 100
par(mfrow = c(2, 4), mgp = c(1.5, 0.5, 0), mar = c(3.25, 3, 1, 0.2) + 0.1) 
plot(B.seq, t0.time[, k], type = "n", ylim = c(min(t1.time[, k]), max(t0.time[, k])), 
     ylab = "time in seconds", xlab = "number of bootstraps", cex.lab = cl, cex.axis = ca)
legend("topleft", text.col = c("blue", "brown", "red", "green"), bty = "n", cex = 0.8,
       legend = c("vectorized multinom. sampl.", 
                  "for loop",
                  "bootstrap R package",
                  "boot R package"))
lines(B.seq, t0.time[, k], col = "red", lwd = 2)
lines(B.seq, t1.time[, k], col = "blue", lwd = 2)
lines(B.seq, t2.time[, k], col = "brown", lwd = 2)
lines(B.seq, t3.time[, k], col = "green", lwd = 2)
mtext("(a)", side = 3, line = 0, at = 9.5e+5, cex = cp)
plot(B.seq, t2.time[, k]/t1.time[, k], col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "time ratio (loop / vect.)", xlab = "number of bootstraps", 
     ylim = c(4, 12))
mtext("(b)", side = 3, line = 0, at = 9.5e+5, cex = cp)
plot(B.seq, t0.time[, k]/t1.time[, k], col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "time ratio (bootstrap / vect.)", xlab = "number of bootstraps", 
     ylim = c(4, 12))
mtext("(c)", side = 3, line = 0, at = 9.5e+5, cex = cp)
plot(B.seq, t3.time[, k]/t1.time[, k], col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "time ratio (boot / vect.)", xlab = "number of bootstraps", 
     ylim = c(4, 12))
mtext("(d)", side = 3, line = 0, at = 9.5e+5, cex = cp)
hist(out1, nclass = nc, probability = TRUE, main = "Vectorized multinomial sampling", xlab = "correlation", 
     cex.lab = cl, cex.axis = ca, xlim = c(0.5, 0.93), cex.main = cm)
mtext("(e)", side = 3, line = -2, at = 0.89, cex = cp)
hist(out2, nclass = nc, probability = TRUE, main = "Resampling (for loop)", xlab = "correlation", 
     cex.lab = cl, cex.axis = ca, xlim = c(0.5, 0.93), cex.main = cm)
mtext("(f)", side = 3, line = -2, at = 0.89, cex = cp)
hist(out0, nclass = nc, probability = TRUE, main = "Resampling (R/bootstrap)", xlab = "correlation", 
     cex.lab = cl, cex.axis = ca, xlim = c(0.5, 0.93), cex.main = cm)
mtext("(g)", side = 3, line = -2, at = 0.89, cex = cp)
hist(out3, nclass = nc, probability = TRUE, main = "Resampling (R/boot)", xlab = "correlation", 
     cex.lab = cl, cex.axis = ca, xlim = c(0.5, 0.93), cex.main = cm)
mtext("(h)", side = 3, line = -2, at = 0.89, cex = cp)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))


#######################################
## American law schools subset example
#######################################

## Same analysis as before, but now focusing on 
## a subset of N = 15 data points

## load the American law schools data subset from the bootstrap R package
##
data(law) 
N <- nrow(law)
x1 <- law$LSAT
x2 <- law$GPA
x <- cbind(x1, x2)

## just check the data-resampling and vectorized implementations 
## on 100,000 bootstraps
##
B <- 100000
system.time(standard1 <- bootstrap(seq(N), B, ThetaCor, x1, x2)$thetastar)
system.time(standard2 <- boot(x, ThetaCor2, B)$t[, 1])
system.time(vectorized <- VectorizedBootstrap(N, B, SampleBootCor, x1, x2))
par(mfrow = c(1, 3))
hist(standard1, nclass = 50, xlim = c(0, 1))
hist(standard2, nclass = 50, xlim = c(0, 1))
hist(vectorized, nclass = 50, xlim = c(0, 1))
par(mfrow = c(1, 1))

## Run the comparison of the vectorized implementation against the 
## data re-sampling implementations in: (i) the bootstrap R package 
## (ii) boot R package; and (iii) on a strait foward for loop 
## implementation.
##
## The implementations are compared over a grid of 30 distinct number 
## of bootstraps values ranging from 1,000 to 1,000,000.
##
B.seq <- round(exp(seq(log(1e+3), log(1e+6), length.out = 30)), 0)
set.seed(123456789)
n.B <- length(B.seq)
t0.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t1.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t2.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t3.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
for(i in seq(n.B)) {
  cat("i = ", i, "\n")
  t0 <- system.time(out0 <- bootstrap(seq(N), B.seq[i], ThetaCor, x1, x2)$thetastar)
  t1 <- system.time(out1 <- VectorizedBootstrap(N, B.seq[i], SampleBootCor, x1, x2))
  t2 <- system.time(out2 <- LoopCorBootstrap(B.seq[i], x1, x2))
  t3 <- system.time(out3 <- boot(x, ThetaCor2, B.seq[i])$t[, 1])
  t0.time[i,] <- t0[1:3]
  t1.time[i,] <- t1[1:3]
  t2.time[i,] <- t2[1:3]
  t3.time[i,] <- t3[1:3]
}

## Generate Figure 2 of the manuscript.
##
cl <- 0.9
ca <- 0.9
cm <- 0.85
cp <- 0.75
k <- 1
nc <- 100
par(mfrow = c(2, 4), mgp = c(1.5, 0.5, 0), mar = c(3.25, 3, 1, 0.2) + 0.1) 
plot(B.seq, t0.time[, k], type = "n", ylim = c(min(t1.time[, k]), max(t0.time[, k])), 
     ylab = "time in seconds", xlab = "number of bootstraps", cex.lab = cl, cex.axis = ca)
legend("topleft", text.col = c("blue", "brown", "red", "green"), bty = "n", cex = 0.8,
       legend = c("vectorized multinom. sampl.", 
                  "for loop",
                  "bootstrap R package",
                  "boot R package"))
lines(B.seq, t0.time[, k], col = "red", lwd = 2)
lines(B.seq, t1.time[, k], col = "blue", lwd = 2)
lines(B.seq, t2.time[, k], col = "brown", lwd = 2)
lines(B.seq, t3.time[, k], col = "green", lwd = 2)
mtext("(a)", side = 3, line = 0, at = 0.5e+5, cex = cp)
plot(B.seq, t2.time[, k]/t1.time[, k], col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "time ratio (loop / vect.)", xlab = "number of bootstraps", 
     ylim = c(5, 90))
mtext("(b)", side = 3, line = 0, at = 0.5e+5, cex = cp)
plot(B.seq, t0.time[, k]/t1.time[, k], col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "time ratio (bootstrap / vect.)", xlab = "number of bootstraps", 
     ylim = c(5, 90))
mtext("(c)", side = 3, line = 0, at = 0.5e+5, cex = cp)
plot(B.seq, t3.time[, k]/t1.time[, k], col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "time ratio (boot / vect.)", xlab = "number of bootstraps", 
     ylim = c(5, 90))
mtext("(d)", side = 3, line = 0, at = 0.5e+5, cex = cp)
hist(out1, nclass = nc, probability = TRUE, main = "Vectorized multinomial sampling", xlab = "correlation", 
     cex.lab = cl, cex.axis = ca, xlim = c(0.1, 1), cex.main = cm)
mtext("(e)", side = 3, line = -2, at = 0.16, cex = cp)
hist(out2, nclass = nc, probability = TRUE, main = "Resampling (for loop)", xlab = "correlation", 
     cex.lab = cl, cex.axis = ca, xlim = c(0.1, 1), cex.main = cm)
mtext("(f)", side = 3, line = -2, at = 0.16, cex = cp)
hist(out0, nclass = nc, probability = TRUE, main = "Resampling (R/bootstrap)", xlab = "correlation", 
     cex.lab = cl, cex.axis = ca, xlim = c(0.1, 1), cex.main = cm)
mtext("(g)", side = 3, line = -2, at = 0.16, cex = cp)
hist(out3, nclass = nc, probability = TRUE, main = "Resampling (R/boot)", xlab = "correlation", 
     cex.lab = cl, cex.axis = ca, xlim = c(0.1, 1), cex.main = cm)
mtext("(h)", side = 3, line = -2, at = 0.16, cex = cp)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))


############################################
## simulated data examples with varying 
## sample sizes
############################################

set.seed(123)
N.max <- 1000
X1 <- rnorm(N.max)
X2 <- X1 + rnorm(N.max)

N.seq <- seq(15, 915, by = 100)
n.N <- length(N.seq)

## 10000 bootstraps
##
B <- 10000
set.seed(987654321)
t0.a <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t1.a <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t2.a <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t3.a <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
for(i in seq(n.N)) {
  cat("i = ", i, "\n")
  x1 <- X1[seq(N.seq[i])]
  x2 <- X2[seq(N.seq[i])]
  x <- cbind(x1, x2)
  t0 <- system.time(out0 <- bootstrap(seq(N.seq[i]), B, ThetaCor, x1, x2)$thetastar)
  t1 <- system.time(out1 <- VectorizedBootstrap(N.seq[i], B, SampleBootCor, x1, x2))
  t2 <- system.time(out2 <- LoopCorBootstrap(B, x1, x2))
  t3 <- system.time(out3 <- boot(x, ThetaCor2, B)$t[, 1])
  t0.a[i,] <- t0[1:3]
  t1.a[i,] <- t1[1:3]
  t2.a[i,] <- t2[1:3]
  t3.a[i,] <- t3[1:3]
}

## 100000 bootstraps
##
B <- 100000
set.seed(123456789)
t0.b <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t1.b <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t2.b <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t3.b <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
for(i in seq(n.N)) {
  cat("i = ", i, "\n")
  x1 <- X1[seq(N.seq[i])]
  x2 <- X2[seq(N.seq[i])]
  x <- cbind(x1, x2)
  t0 <- system.time(out0 <- bootstrap(seq(N.seq[i]), B, ThetaCor, x1, x2)$thetastar)
  t1 <- system.time(out1 <- VectorizedBootstrap(N.seq[i], B, SampleBootCor, x1, x2))
  t2 <- system.time(out2 <- LoopCorBootstrap(B, x1, x2))
  t3 <- system.time(out3 <- boot(x, ThetaCor2, B)$t[, 1])
  t0.b[i,] <- t0[1:3]
  t1.b[i,] <- t1[1:3]
  t2.b[i,] <- t2[1:3]
  t3.b[i,] <- t3[1:3]
}

## 1000000 bootstraps
##
B <- 1000000
set.seed(123456789)
t0.c <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t1.c <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t2.c <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t3.c <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
for(i in seq(n.N)) {
  cat("i = ", i, "\n")
  x1 <- X1[seq(N.seq[i])]
  x2 <- X2[seq(N.seq[i])]
  x <- cbind(x1, x2)
  t0 <- system.time(out0 <- bootstrap(seq(N.seq[i]), B, ThetaCor, x1, x2)$thetastar)
  t1 <- system.time(out1 <- VectorizedBootstrap(N.seq[i], B, SampleBootCor, x1, x2))
  t2 <- system.time(out2 <- LoopCorBootstrap(B, x1, x2))
  t3 <- system.time(out3 <- boot(x, ThetaCor2, B)$t[, 1])
  t0.c[i,] <- t0[1:3]
  t1.c[i,] <- t1[1:3]
  t2.c[i,] <- t2[1:3]
  t3.c[i,] <- t3[1:3]
}

## Generate Figure 3 of the manuscript.
##
cl <- 0.8
ca <- 0.8
cm <- 0.9
cp <- 0.8
k <- 1
par(mfrow = c(3, 4), mgp = c(1.5, 0.5, 0), mar = c(2.5, 3, 1, 0.2) + 0.1)
##
plot(N.seq, t0.a[, k], type = "n", ylim = c(min(t1.a[, k]), max(t0.a[, k])), 
     ylab = "time in seconds", xlab = "sample size", cex.lab = cl, cex.axis = ca,
     main = "B = 10,000", las = 2, cex.main = cm)
legend("bottomright", text.col = c("blue", "brown", "red", "green"), bty = "n", cex = 0.8,
       legend = c("vectorized", "for loop", "R/bootstrap", "R/boot"))
lines(N.seq, t0.a[, k], col = "red", lwd = 2)
lines(N.seq, t1.a[, k], col = "blue", lwd = 2)
lines(N.seq, t2.a[, k], col = "brown", lwd = 2)
lines(N.seq, t3.a[, k], col = "green", lwd = 2)
mtext("(a)", side = 3, line = 0, at = 900, cex = cp)
plot(N.seq, log(t2.a[, k]/t1.a[, k]), col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "log time ratio (loop / vect.)", xlab = "sample size", cex.main = cm, 
     ylim = c(-0.2, 4.2), main = "B = 10,000", type = "b", lwd = 2, las = 2)
abline(h = 0)
mtext("(b)", side = 3, line = 0, at = 900, cex = cp)
plot(N.seq, log(t0.a[, k]/t1.a[, k]), col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "log time ratio (bootstrap / vect.)", xlab = "sample size", cex.main = cm, 
     ylim = c(-0.2, 4.2), main = "B = 10,000", type = "b", lwd = 2, las = 2)
abline(h = 0)
mtext("(c)", side = 3, line = 0, at = 900, cex = cp)
plot(N.seq, log(t3.a[, k]/t1.a[, k]), col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "log time ratio (boot / vect.)", xlab = "sample size", cex.main = cm, 
     ylim = c(-0.2, 4.2), main = "B = 10,000", type = "b", lwd = 2, las = 2)
abline(h = 0)
mtext("(d)", side = 3, line = 0, at = 900, cex = cp)
##
plot(N.seq, t0.b[, k], type = "n", ylim = c(min(t1.b[, k]), max(t0.b[, k])), 
     ylab = "time in seconds", xlab = "sample size", cex.lab = cl, cex.axis = ca,
     main = "B = 100,000", las = 2, cex.main = cm)
legend("bottomright", text.col = c("blue", "brown", "red", "green"), bty = "n", cex = 0.8,
       legend = c("vectorized", "for loop", "R/bootstrap", "R/boot"))
lines(N.seq, t0.b[, k], col = "red", lwd = 2)
lines(N.seq, t1.b[, k], col = "blue", lwd = 2)
lines(N.seq, t2.b[, k], col = "brown", lwd = 2)
lines(N.seq, t3.b[, k], col = "green", lwd = 2)
mtext("(e)", side = 3, line = 0, at = 900, cex = cp)
plot(N.seq, log(t2.b[, k]/t1.b[, k]), col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "log time ratio (loop / vect.)", xlab = "sample size", cex.main = cm, 
     ylim = c(-0.2, 4.2), main = "B = 100,000", type = "b", lwd = 2, las = 2)
abline(h = 0)
mtext("(f)", side = 3, line = 0, at = 900, cex = cp)
plot(N.seq, log(t0.b[, k]/t1.b[, k]), col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "log time ratio (bootstrap / vect.)", xlab = "sample size", cex.main = cm,
     ylim = c(-0.2, 4.2), main = "B = 100,000", type = "b", lwd = 2, las = 2)
abline(h = 0)
mtext("(g)", side = 3, line = 0, at = 900, cex = cp)
plot(N.seq, log(t3.b[, k]/t1.b[, k]), col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "log time ratio (boot / vect.)", xlab = "sample size", cex.main = cm, 
     ylim = c(-0.2, 4.2), main = "B = 100,000", type = "b", lwd = 2, las = 2)
abline(h = 0)
mtext("(h)", side = 3, line = 0, at = 900, cex = cp)
##
plot(N.seq, t0.c[, k], type = "n", ylim = c(min(t1.c[, k]), max(t0.c[, k])), 
     ylab = "time in seconds", xlab = "sample size", cex.lab = cl, cex.axis = ca,
     main = "B = 1,000,000", las = 2, cex.main = cm)
legend("bottomright", text.col = c("blue", "brown", "red", "green"), bty = "n", cex = 0.8,
       legend = c("vectorized", "for loop", "R/bootstrap", "R/boot"))
lines(N.seq, t0.c[, k], col = "red", lwd = 2)
lines(N.seq, t1.c[, k], col = "blue", lwd = 2)
lines(N.seq, t2.c[, k], col = "brown", lwd = 2)
lines(N.seq, t3.c[, k], col = "green", lwd = 2)
mtext("(i)", side = 3, line = 0, at = 900, cex = cp)
plot(N.seq, log(t2.c[, k]/t1.c[, k]), col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "log time ratio (loop / vect.)", xlab = "sample size", cex.main = cm, 
     ylim = c(-0.2, 4.2), main = "B = 1,000,000", type = "b", lwd = 2, las = 2)
abline(h = 0)
mtext("(j)", side = 3, line = 0, at = 900, cex = cp)
plot(N.seq, log(t0.c[, k]/t1.c[, k]), col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "log time ratio (bootstrap / vect.)", xlab = "sample size", cex.main = cm, 
     ylim = c(-0.2, 4.2), main = "B = 1,000,000", type = "b", lwd = 2, las = 2)
abline(h = 0)
mtext("(k)", side = 3, line = 0, at = 900, cex = cp)
plot(N.seq, log(t3.c[, k]/t1.c[, k]), col = "black", cex.lab = cl, cex.axis = ca,
     ylab = "log time ratio (boot / vect.)", xlab = "sample size", cex.main = cm, 
     ylim = c(-0.2, 4.2), main = "B = 1,000,000", type = "b", lwd = 2, las = 2)
abline(h = 0)
mtext("(l)", side = 3, line = 0, at = 900, cex = cp)
par(mfrow = c(1, 1), mgp = c(3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)


## Generate the results for Fig. 4
##
set.seed(123)
N.max <- 1000
X1 <- rnorm(N.max)
X2 <- X1 + rnorm(N.max)

N.seq <- seq(15, 915, by = 100)
n.N <- length(N.seq)

set.seed(123456789)
out <- ReplicateVecTiming(K = 10, N.seq, x1, x2)

t1.W <- out$at1.W
t2.W <- out$at2.W
t3.W <- out$at3.W

t1.M <- out$at1.M
t2.M <- out$at2.M
t3.M <- out$at3.M

## Generate Fig. 4 of the manuscript
##
k <- 1
cp <- 0.8
par(mfrow = c(1, 4), mgp = c(1.5, 0.5, 0), mar = c(3.25, 3, 1, 0.3) + 0.1) ##c(bottom, left, top, right)
plot(N.seq, t1.W[, k], type = "l", lwd = 2, col = "blue", ylab = "time in seconds",
     xlab = "sample size", main = "B = 10,000")
lines(N.seq, t1.M[, k], lwd = 2, col = "red")
mtext("(a)", side = 3, line = 0, at = 900, cex = cp)
legend("topleft", legend = c("multinomial sampler", "matrix operations"), 
       text.col = c("blue", "red"), bty = "n", cex = 0.9)
plot(N.seq, t2.W[, k], type = "l", lwd = 2, col = "blue", ylab = "time in seconds",
     xlab = "sample size", main = "B = 100,000")
lines(N.seq, t2.M[, k], lwd = 2, col = "red")
mtext("(b)", side = 3, line = 0, at = 900, cex = cp)
legend("topleft", legend = c("multinomial sampler", "matrix operations"), 
       text.col = c("blue", "red"), bty = "n", cex = 0.9)
plot(N.seq, t3.W[, k], type = "l", lwd = 2, col = "blue", ylab = "time in seconds",
     xlab = "sample size", main = "B = 1,000,000")
lines(N.seq, t3.M[, k], lwd = 2, col = "red")
mtext("(c)", side = 3, line = 0, at = 900, cex = cp)
legend("topleft", legend = c("multinomial sampler", "matrix operations"), 
       text.col = c("blue", "red"), bty = "n", cex = 0.9)
plot(N.seq, t3.W[, k]/t3.M[, k], lwd = 2, col = "black", ylab = "time ratio in seconds",
     xlab = "sample size", ylim = c(0, 20), type = "l")
lines(N.seq, t2.W[, k]/t2.M[, k], lwd = 2, col = "darkgreen")
lines(N.seq, t1.W[, k]/t1.M[, k], lwd = 2, col = "orange3")
legend("bottomright", legend = c("B = 1,000,000", "B = 100,000", "B = 10,000"), 
       text.col = c("black", "darkgreen", "orange"), bty = "n")
mtext("(d)", side = 3, line = 0, at = 900, cex = cp)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))
