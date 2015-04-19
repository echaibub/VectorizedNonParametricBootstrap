#######################################################
## Script generating the results in Supplementary file
## S2 Text.
#######################################################

library(devtools)
source_url("https://raw.githubusercontent.com/echaibub/VectorizedNonParametricBootstrap/master/vectorized_bootstrap_functions.R")
source_url("https://raw.githubusercontent.com/echaibub/VectorizedNonParametricBootstrap/master/parallel_cor_boot_functions.R")

###################################################
###################################################
###################################################

## we run these scripts in the Windows and the Xubuntu machines

data(law82) 
N <- nrow(law82)
x1 <- law82$LSAT
x2 <- law82$GPA
x <- cbind(x1, x2)

B.seq <- round(exp(seq(log(1e+3), log(1e+6), length.out = 30)), 0)
set.seed(123456789)
n.B <- length(B.seq)
t1.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t2.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t3.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t4.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t5.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t6.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
for(i in seq(n.B)) {
  cat("law82 = ", i, "\n")
  t1 <- system.time(out1 <- LoopCorBootstrap(B.seq[i], x1, x2))
  t2 <- system.time(out2 <- ParallelLoopCorBootstrap(B.seq[i], x1, x2))
  t3 <- system.time(out3 <- ParallelLoopCorBootstrap2(B.seq[i], x1, x2))
  t4 <- system.time(out4 <- VectorizedBootstrap(N, B.seq[i], SampleBootCor, x1, x2))
  t5 <- system.time(out5 <- ParallelVectorizedCorBootstrap(B.seq[i], N, x1, x2))
  t6 <- system.time(out6 <- ParallelVectorizedCorBootstrap2(B.seq[i], N, x1, x2))
  t1.time[i,] <- t1[1:3]
  t2.time[i,] <- t2[1:3]
  t3.time[i,] <- t3[1:3]
  t4.time[i,] <- t4[1:3]
  t5.time[i,] <- t5[1:3]
  t6.time[i,] <- t6[1:3]
}
#save(t1.time, t2.time, t3.time, t4.time, t5.time, t6.time,
#     file = "law82_parallel_time_comparison.RData", compress = TRUE)

#################################################
#################################################
#################################################

data(law) 
N <- nrow(law)
x1 <- law$LSAT
x2 <- law$GPA
x <- cbind(x1, x2)

B.seq <- round(exp(seq(log(1e+3), log(1e+6), length.out = 30)), 0)
set.seed(123456789)
n.B <- length(B.seq)
t1.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t2.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t3.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t4.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t5.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
t6.time <- matrix(NA, n.B, 3, dimnames = list(as.character(B.seq), c("user", "system", "elapsed")))
for(i in seq(n.B)) {
  cat("law = ", i, "\n")
  t1 <- system.time(out1 <- LoopCorBootstrap(B.seq[i], x1, x2))
  t2 <- system.time(out2 <- ParallelLoopCorBootstrap(B.seq[i], x1, x2))
  t3 <- system.time(out3 <- ParallelLoopCorBootstrap2(B.seq[i], x1, x2))
  t4 <- system.time(out4 <- VectorizedBootstrap(N, B.seq[i], SampleBootCor, x1, x2))
  t5 <- system.time(out5 <- ParallelVectorizedCorBootstrap(B.seq[i], N, x1, x2))
  t6 <- system.time(out6 <- ParallelVectorizedCorBootstrap2(B.seq[i], N, x1, x2))
  t1.time[i,] <- t1[1:3]
  t2.time[i,] <- t2[1:3]
  t3.time[i,] <- t3[1:3]
  t4.time[i,] <- t4[1:3]
  t5.time[i,] <- t5[1:3]
  t6.time[i,] <- t6[1:3]
}
#save(t1.time, t2.time, t3.time, t4.time, t5.time, t6.time,
#     file = "law_parallel_time_comparison.RData", compress = TRUE)


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
t1.a <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t2.a <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t3.a <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t4.a <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t5.a <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t6.a <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
for(i in seq(n.N)) {
  cat("B = 10000 ", i, "\n")
  x1 <- X1[seq(N.seq[i])]
  x2 <- X2[seq(N.seq[i])]
  x <- cbind(x1, x2)
  t1 <- system.time(out1 <- LoopCorBootstrap(B, x1, x2))
  t2 <- system.time(out2 <- ParallelLoopCorBootstrap(B, x1, x2))
  t3 <- system.time(out3 <- ParallelLoopCorBootstrap2(B, x1, x2))
  t4 <- system.time(out4 <- VectorizedBootstrap(N.seq[i], B, SampleBootCor, x1, x2))
  t5 <- system.time(out5 <- ParallelVectorizedCorBootstrap(B, N.seq[i], x1, x2)) 
  t6 <- system.time(out6 <- ParallelVectorizedCorBootstrap2(B, N.seq[i], x1, x2))  
  t1.a[i,] <- t1[1:3]
  t2.a[i,] <- t2[1:3]
  t3.a[i,] <- t3[1:3]
  t4.a[i,] <- t4[1:3]
  t5.a[i,] <- t5[1:3]
  t6.a[i,] <- t6[1:3]
}
#save(t1.a, t2.a, t3.a, t4.a, t5.a, t6.a, file = "parallel_boot_times_10000B.RData", compress = TRUE)



## 100000 bootstraps
##
B <- 100000
set.seed(987654321)
t1.b <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t2.b <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t3.b <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t4.b <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t5.b <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t6.b <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
for(i in seq(n.N)) {
  cat("B = 100000 ", i, "\n")
  x1 <- X1[seq(N.seq[i])]
  x2 <- X2[seq(N.seq[i])]
  x <- cbind(x1, x2)
  t1 <- system.time(out1 <- LoopCorBootstrap(B, x1, x2))
  t2 <- system.time(out2 <- ParallelLoopCorBootstrap(B, x1, x2))
  t3 <- system.time(out3 <- ParallelLoopCorBootstrap2(B, x1, x2))
  t4 <- system.time(out4 <- VectorizedBootstrap(N.seq[i], B, SampleBootCor, x1, x2))
  t5 <- system.time(out5 <- ParallelVectorizedCorBootstrap(B, N.seq[i], x1, x2)) 
  t6 <- system.time(out6 <- ParallelVectorizedCorBootstrap2(B, N.seq[i], x1, x2))  
  t1.b[i,] <- t1[1:3]
  t2.b[i,] <- t2[1:3]
  t3.b[i,] <- t3[1:3]
  t4.b[i,] <- t4[1:3]
  t5.b[i,] <- t5[1:3]
  t6.b[i,] <- t6[1:3]
}
#save(t1.b, t2.b, t3.b, t4.b, t5.b, t6.b, file = "parallel_boot_times_100000B.RData", compress = TRUE)


## 1000000 bootstraps
##
B <- 1000000
set.seed(987654321)
t1.c <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t2.c <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t3.c <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t4.c <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t5.c <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
t6.c <- matrix(NA, n.N, 3, dimnames = list(as.character(N.seq), c("user", "system", "elapsed")))
for(i in seq(n.N)) {
  cat("B = 1000000", i, "\n")
  x1 <- X1[seq(N.seq[i])]
  x2 <- X2[seq(N.seq[i])]
  x <- cbind(x1, x2)
  t1 <- system.time(out1 <- LoopCorBootstrap(B, x1, x2))
  t2 <- system.time(out2 <- ParallelLoopCorBootstrap(B, x1, x2))
  t3 <- system.time(out3 <- ParallelLoopCorBootstrap2(B, x1, x2))
  t4 <- system.time(out4 <- VectorizedBootstrap(N.seq[i], B, SampleBootCor, x1, x2))
  t5 <- system.time(out5 <- ParallelVectorizedCorBootstrap(B, N.seq[i], x1, x2))
  t6 <- system.time(out6 <- ParallelVectorizedCorBootstrap2(B, N.seq[i], x1, x2))  
  t1.c[i,] <- t1[1:3]
  t2.c[i,] <- t2[1:3]
  t3.c[i,] <- t3[1:3]
  t4.c[i,] <- t4[1:3]
  t5.c[i,] <- t5[1:3]
  t6.c[i,] <- t6[1:3]
}
#save(t1.c, t2.c, t3.c, t4.c, t5.c, t6.c, file = "parallel_boot_times_1000000B.RData", compress = TRUE)




##########################################################################
## generates Figures 1, 2, 3, 4, and 5 in Text S2
##########################################################################

###########################################
## law82 results
###########################################

B.seq <- round(exp(seq(log(1e+3), log(1e+6), length.out = 30)), 0)

## law82

cl <- 1.3
ca <- 1.3
cp <- 1.3
cle <- 0.9
k <- 3

load("law82_parallel_time_comparison_windows.RData")
tmin1 <- min(c(t1.time[, k], t2.time[, k], t3.time[, k], t4.time[, k],  t5.time[, k], t6.time[, k]))
tmax1 <- max(c(t1.time[, k], t2.time[, k], t3.time[, k], t4.time[, k],  t5.time[, k], t6.time[, k]))

load("law82_parallel_time_comparison_ubuntu.RData")
tmin2 <- min(c(t1.time[, k], t2.time[, k], t3.time[, k], t4.time[, k],  t5.time[, k], t6.time[, k]))
tmax2 <- max(c(t1.time[, k], t2.time[, k], t3.time[, k], t4.time[, k],  t5.time[, k], t6.time[, k]))

tmin <- min(c(tmin1, tmin2))
tmax <- max(c(tmax1, tmax2))


par(mfrow = c(1, 2), mar = c(4, 4, 3, 2) + 0.1)
###
load("law82_parallel_time_comparison_windows.RData")
plot(B.seq, t1.time[, k], type = "n", ylim = c(tmin, tmax), main = "Windows,  N = 82",
     ylab = "time in seconds", xlab = "number of bootstraps", cex.lab = cl, cex.axis = ca)
legend("topleft", text.col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       col = c("brown", "brown", "brown", "blue", "blue", "blue"),
       lty = c(1, 2, 3, 1, 2, 3), bty = "n", cex = cle,
       legend = c("for loop", 
                  "parLapply - for loop",
                  "mclapply - for loop",
                  "vectorized",
                  "parLapply - vectorized",
                  "mclapply - vectorized"))
lines(B.seq, t1.time[, k], col = "brown", lwd = 2)
lines(B.seq, t2.time[, k], col = "brown", lwd = 2, lty = 2)
lines(B.seq, t3.time[, k], col = "brown", lwd = 2, lty = 3)
lines(B.seq, t4.time[, k], col = "blue", lwd = 2)
lines(B.seq, t5.time[, k], col = "blue", lwd = 2, lty = 2)
lines(B.seq, t6.time[, k], col = "blue", lwd = 2, lty = 3)
mtext("(a)", side = 3, line = 0.5, at = 9.5e+5, cex = cp)
###
load("law82_parallel_time_comparison_ubuntu.RData")
plot(B.seq, t1.time[, k], type = "n", ylim = c(tmin, tmax), main = "Xubuntu,  N = 82",
     ylab = "time in seconds", xlab = "number of bootstraps", cex.lab = cl, cex.axis = ca)
legend("topleft", text.col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       col = c("brown", "brown", "brown", "blue", "blue", "blue"),
       lty = c(1, 2, 3, 1, 2, 3), bty = "n", cex = cle,
       legend = c("for loop", 
                  "parLapply - for loop",
                  "mclapply - for loop",
                  "vectorized",
                  "parLapply - vectorized",
                  "mclapply - vectorized"))
lines(B.seq, t1.time[, k], col = "brown", lwd = 2)
lines(B.seq, t2.time[, k], col = "brown", lwd = 2, lty = 2)
lines(B.seq, t3.time[, k], col = "brown", lwd = 2, lty = 3)
lines(B.seq, t4.time[, k], col = "blue", lwd = 2)
lines(B.seq, t5.time[, k], col = "blue", lwd = 2, lty = 2)
lines(B.seq, t6.time[, k], col = "blue", lwd = 2, lty = 3)
mtext("(b)", side = 3, line = 0.5, at = 9.5e+5, cex = cp)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)



###########################################
## law results
###########################################

outpath <- "C:/Users/echaibub/Documents/FastBootstrap/Plos revision/"

B.seq <- round(exp(seq(log(1e+3), log(1e+6), length.out = 30)), 0)

## law

cl <- 1.3
ca <- 1.3
cp <- 1.3
cle <- 0.9
k <- 3

load("law_parallel_time_comparison_windows.RData")
tmin1 <- min(c(t1.time[, k], t2.time[, k], t3.time[, k], t4.time[, k],  t5.time[, k], t6.time[, k]))
tmax1 <- max(c(t1.time[, k], t2.time[, k], t3.time[, k], t4.time[, k],  t5.time[, k], t6.time[, k]))

load("law_parallel_time_comparison_ubuntu.RData")
tmin2 <- min(c(t1.time[, k], t2.time[, k], t3.time[, k], t4.time[, k],  t5.time[, k], t6.time[, k]))
tmax2 <- max(c(t1.time[, k], t2.time[, k], t3.time[, k], t4.time[, k],  t5.time[, k], t6.time[, k]))

tmin <- min(c(tmin1, tmin2))
tmax <- max(c(tmax1, tmax2))


par(mfrow = c(1, 2), mar = c(4, 4, 3, 2) + 0.1)
###
load("law_parallel_time_comparison_windows.RData")
plot(B.seq, t1.time[, k], type = "n", ylim = c(tmin, tmax), main = "Windows,  N = 15",
     ylab = "time in seconds", xlab = "number of bootstraps", cex.lab = cl, cex.axis = ca)
legend("topleft", text.col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       col = c("brown", "brown", "brown", "blue", "blue", "blue"),
       lty = c(1, 2, 3, 1, 2, 3), bty = "n", cex = cle,
       legend = c("for loop", 
                  "parLapply - for loop",
                  "mclapply - for loop",
                  "vectorized",
                  "parLapply - vectorized",
                  "mclapply - vectorized"))
lines(B.seq, t1.time[, k], col = "brown", lwd = 2)
lines(B.seq, t2.time[, k], col = "brown", lwd = 2, lty = 2)
lines(B.seq, t3.time[, k], col = "brown", lwd = 2, lty = 3)
lines(B.seq, t4.time[, k], col = "blue", lwd = 2)
lines(B.seq, t5.time[, k], col = "blue", lwd = 2, lty = 2)
lines(B.seq, t6.time[, k], col = "blue", lwd = 2, lty = 3)
mtext("(a)", side = 3, line = 0.5, at = 9.5e+5, cex = cp)
###
load("law_parallel_time_comparison_ubuntu.RData")
plot(B.seq, t1.time[, k], type = "n", ylim = c(tmin, tmax), main = "Xubuntu,  N = 15",
     ylab = "time in seconds", xlab = "number of bootstraps", cex.lab = cl, cex.axis = ca)
legend("topleft", text.col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       col = c("brown", "brown", "brown", "blue", "blue", "blue"),
       lty = c(1, 2, 3, 1, 2, 3), bty = "n", cex = cle,
       legend = c("for loop", 
                  "parLapply - for loop",
                  "mclapply - for loop",
                  "vectorized",
                  "parLapply - vectorized",
                  "mclapply - vectorized"))
lines(B.seq, t1.time[, k], col = "brown", lwd = 2)
lines(B.seq, t2.time[, k], col = "brown", lwd = 2, lty = 2)
lines(B.seq, t3.time[, k], col = "brown", lwd = 2, lty = 3)
lines(B.seq, t4.time[, k], col = "blue", lwd = 2)
lines(B.seq, t5.time[, k], col = "blue", lwd = 2, lty = 2)
lines(B.seq, t6.time[, k], col = "blue", lwd = 2, lty = 3)
mtext("(b)", side = 3, line = 0.5, at = 9.5e+5, cex = cp)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)



###########################################
## simulated results
###########################################

N.seq <- seq(15, 915, by = 100)

cl <- 1.3
ca <- 1.3
cp <- 1.3
cle <- 0.9
k <- 3

## B = 10000, seq.N

load("parallel_boot_times_10000B_windows.RData")
amin1 <- min(c(t1.a[, k], t2.a[, k], t3.a[, k], t4.a[, k], t5.a[, k], t6.a[, k]))
amax1 <- max(c(t1.a[, k], t2.a[, k], t3.a[, k], t4.a[, k], t5.a[, k], t6.a[, k]))

load("parallel_boot_times_10000B_ubuntu.RData")
amin2 <- min(c(t1.a[, k], t2.a[, k], t3.a[, k], t4.a[, k], t5.a[, k], t6.a[, k]))
amax2 <- max(c(t1.a[, k], t2.a[, k], t3.a[, k], t4.a[, k], t5.a[, k], t6.a[, k]))

amin <- min(c(amin1, amin2))
amax <- max(c(amax1, amax2))


par(mfrow = c(1, 2), mar = c(4, 4, 3, 2) + 0.1)
###
load("parallel_boot_times_10000B_windows.RData")
plot(N.seq, t1.a[, k], type = "n", ylim = c(amin, amax), 
     ylab = "time in seconds", xlab = "sample size", cex.lab = cl, cex.axis = ca,
     main = "Windows,  B = 10,000", las = 1)
legend("topleft", text.col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       lty = c(1, 2, 3, 1, 2, 3), bty = "n", cex = cle,
       legend = c("for loop", "parLapply - for loop", "mclapply- for loop", 
                  "vectorized", "parLapply - vectorized", "mclapply - vectorized"))
lines(N.seq, t1.a[, k], col = "brown", lwd = 2)
lines(N.seq, t2.a[, k], col = "brown", lwd = 2, lty = 2)
lines(N.seq, t3.a[, k], col = "brown", lwd = 2, lty = 3)
lines(N.seq, t4.a[, k], col = "blue", lwd = 2)
lines(N.seq, t5.a[, k], col = "blue", lwd = 2, lty = 2)
lines(N.seq, t6.a[, k], col = "blue", lwd = 2, lty = 3)
mtext("(a)", side = 3, line = 0.5, at = 900, cex = cp)
###
load("parallel_boot_times_10000B_ubuntu.RData")
plot(N.seq, t1.a[, k], type = "n", ylim = c(amin, amax), 
     ylab = "time in seconds", xlab = "sample size", cex.lab = cl, cex.axis = ca,
     main = "Xubuntu,  B = 10,000", las = 1)
legend("bottomright", text.col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       lty = c(1, 2, 3, 1, 2, 3), bty = "n", cex = 0.8,
       legend = c("for loop", "parLapply - for loop", "mclapply- for loop", 
                  "vectorized", "parLapply - vectorized", "mclapply - vectorized"))
lines(N.seq, t1.a[, k], col = "brown", lwd = 2)
lines(N.seq, t2.a[, k], col = "brown", lwd = 2, lty = 2)
lines(N.seq, t3.a[, k], col = "brown", lwd = 2, lty = 3)
lines(N.seq, t4.a[, k], col = "blue", lwd = 2)
lines(N.seq, t5.a[, k], col = "blue", lwd = 2, lty = 2)
lines(N.seq, t6.a[, k], col = "blue", lwd = 2, lty = 3)
mtext("(b)", side = 3, line = 0.5, at = 900, cex = cp)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)




## B = 100000, seq.N

load("parallel_boot_times_100000B_windows.RData")
bmin1 <- min(c(t1.b[, k], t2.b[, k], t3.b[, k], t4.b[, k], t5.b[, k], t6.b[, k]))
bmax1 <- max(c(t1.b[, k], t2.b[, k], t3.b[, k], t4.b[, k], t5.b[, k], t6.b[, k]))

load("parallel_boot_times_100000B_ubuntu.RData")
bmin2 <- min(c(t1.b[, k], t2.b[, k], t3.b[, k], t4.b[, k], t5.b[, k], t6.b[, k]))
bmax2 <- max(c(t1.b[, k], t2.b[, k], t3.b[, k], t4.b[, k], t5.b[, k], t6.b[, k]))

bmin <- min(c(bmin1, bmin2))
bmax <- max(c(bmax1, bmax2))


par(mfrow = c(1, 2), mar = c(4, 4, 3, 2) + 0.1)
###
load("parallel_boot_times_100000B_windows.RData")
plot(N.seq, t1.b[, k], type = "n", ylim = c(bmin, bmax), 
     ylab = "time in seconds", xlab = "sample size", cex.lab = cl, cex.axis = ca,
     main = "Windows,  B = 100,000", las = 1)
legend("topleft", text.col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       lty = c(1, 2, 3, 1, 2, 3), bty = "n", cex = cle,
       legend = c("for loop", "parLapply - for loop", "mclapply- for loop", 
                  "vectorized", "parLapply - vectorized", "mclapply - vectorized"))
lines(N.seq, t1.b[, k], col = "brown", lwd = 2)
lines(N.seq, t2.b[, k], col = "brown", lwd = 2, lty = 2)
lines(N.seq, t3.b[, k], col = "brown", lwd = 2, lty = 3)
lines(N.seq, t4.b[, k], col = "blue", lwd = 2)
lines(N.seq, t5.b[, k], col = "blue", lwd = 2, lty = 2)
lines(N.seq, t6.b[, k], col = "blue", lwd = 2, lty = 3)
mtext("(a)", side = 3, line = 0.5, at = 900, cex = cp)
###
load("parallel_boot_times_100000B_ubuntu.RData")
plot(N.seq, t1.b[, k], type = "n", ylim = c(bmin, bmax), 
     ylab = "time in seconds", xlab = "sample size", cex.lab = cl, cex.axis = ca,
     main = "Xubuntu,  B = 100,000", las = 1)
legend("topleft", text.col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       lty = c(1, 2, 3, 1, 2, 3), bty = "n", cex = cle,
       legend = c("for loop", "parLapply - for loop", "mclapply- for loop", 
                  "vectorized", "parLapply - vectorized", "mclapply - vectorized"))
lines(N.seq, t1.b[, k], col = "brown", lwd = 2)
lines(N.seq, t2.b[, k], col = "brown", lwd = 2, lty = 2)
lines(N.seq, t3.b[, k], col = "brown", lwd = 2, lty = 3)
lines(N.seq, t4.b[, k], col = "blue", lwd = 2)
lines(N.seq, t5.b[, k], col = "blue", lwd = 2, lty = 2)
lines(N.seq, t6.b[, k], col = "blue", lwd = 2, lty = 3)
mtext("(b)", side = 3, line = 0.5, at = 900, cex = cp)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)



## B = 1000000, seq.N

load("parallel_boot_times_1000000B_windows.RData")
cmin1 <- min(c(t1.c[, k], t2.c[, k], t3.c[, k], t4.c[, k], t5.c[, k], t6.c[, k]))
cmax1 <- max(c(t1.c[, k], t2.c[, k], t3.c[, k], t4.c[, k], t5.c[, k], t6.c[, k]))

load("parallel_boot_times_1000000B_ubuntu.RData")
cmin2 <- min(c(t1.c[, k], t2.c[, k], t3.c[, k], t4.c[, k], t5.c[, k], t6.c[, k]))
cmax2 <- max(c(t1.c[, k], t2.c[, k], t3.c[, k], t4.c[, k], t5.c[, k], t6.c[, k]))

cmin <- min(c(cmin1, cmin2))
cmax <- max(c(cmax1, cmax2))


par(mfrow = c(1, 2), mar = c(4, 4, 3, 2) + 0.1)
###
load("parallel_boot_times_1000000B_windows.RData")
plot(N.seq, t1.c[, k], type = "n", ylim = c(cmin, cmax), 
     ylab = "time in seconds", xlab = "sample size", cex.lab = cl, cex.axis = ca,
     main = "Windows,  B = 1,000,000", las = 1)
legend("topleft", text.col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       lty = c(1, 2, 3, 1, 2, 3), bty = "n", cex = cle,
       legend = c("for loop", "parLapply - for loop", "mclapply- for loop", 
                  "vectorized", "parLapply - vectorized", "mclapply - vectorized"))
lines(N.seq, t1.c[, k], col = "brown", lwd = 2)
lines(N.seq, t2.c[, k], col = "brown", lwd = 2, lty = 2)
lines(N.seq, t3.c[, k], col = "brown", lwd = 2, lty = 3)
lines(N.seq, t4.c[, k], col = "blue", lwd = 2)
lines(N.seq, t5.c[, k], col = "blue", lwd = 2, lty = 2)
lines(N.seq, t6.c[, k], col = "blue", lwd = 2, lty = 3)
mtext("(a)", side = 3, line = 0.5, at = 900, cex = cp)
###
load("parallel_boot_times_1000000B_ubuntu.RData")
plot(N.seq, t1.c[, k], type = "n", ylim = c(cmin, cmax), 
     ylab = "time in seconds", xlab = "sample size", cex.lab = cl, cex.axis = ca,
     main = "Xubuntu,  B = 1,000,000", las = 1)
legend("topleft", text.col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       col = c("brown", "brown", "brown", "blue", "blue", "blue"), 
       lty = c(1, 2, 3, 1, 2, 3), bty = "n", cex = cle,
       legend = c("for loop", "parLapply - for loop", "mclapply- for loop", 
                  "vectorized", "parLapply - vectorized", "mclapply - vectorized"))
lines(N.seq, t1.c[, k], col = "brown", lwd = 2)
lines(N.seq, t2.c[, k], col = "brown", lwd = 2, lty = 2)
lines(N.seq, t3.c[, k], col = "brown", lwd = 2, lty = 3)
lines(N.seq, t4.c[, k], col = "blue", lwd = 2)
lines(N.seq, t5.c[, k], col = "blue", lwd = 2, lty = 2)
lines(N.seq, t6.c[, k], col = "blue", lwd = 2, lty = 3)
mtext("(b)", side = 3, line = 0.5, at = 900, cex = cp)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

