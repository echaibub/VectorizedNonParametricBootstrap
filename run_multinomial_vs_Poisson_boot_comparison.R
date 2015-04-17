
###########################################################
## Run the comparison between the vectorized multinomial 
## bootstrap and the vectorized implementation of the
## approximate bootstrap based on Poisson frequencies.
###########################################################

library(devtools)
source_url("https://raw.githubusercontent.com/echaibub/VectorizedNonParametricBootstrap/master/vectorized_bootstrap_functions.R")
source_url("https://raw.githubusercontent.com/echaibub/VectorizedNonParametricBootstrap/master/vectorized_Poisson_frequencies_functions.R")

############################################################
## run the comparisons on the law82 data (N = 82)
############################################################

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

for(i in seq(n.B)) {
  cat("law82 = ", i, "\n")
  t1 <- system.time(out1 <- VectorizedBootstrap(N, B.seq[i], SampleBootCor, x1, x2))
  t2 <- system.time(out2 <- VectorizedPoissonBootstrap(N, B.seq[i], SampleBootCor2, x1, x2))
  t1.time[i,] <- t1[1:3]
  t2.time[i,] <- t2[1:3]
}

#save(t1.time, t2.time, out1, out2,
#     file = "law82_poisson_comparison.RData", compress = TRUE)



############################################################
## run the comparisons on the law data (subset with N = 15)
############################################################

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

for(i in seq(n.B)) {
  cat("law82 = ", i, "\n")
  t1 <- system.time(out1 <- VectorizedBootstrap(N, B.seq[i], SampleBootCor, x1, x2))
  t2 <- system.time(out2 <- VectorizedPoissonBootstrap(N, B.seq[i], SampleBootCor2, x1, x2))
  t1.time[i,] <- t1[1:3]
  t2.time[i,] <- t2[1:3]
}

#save(t1.time, t2.time, out1, out2,
#     file = "law_poisson_comparison.RData", compress = TRUE)


##########################################################
## Generate Figures 1 and 2 of Supplementary file S3 Text
##########################################################

## Figure 1

cl <- 1.3
ca <- 1.3
cp <- 1.3
cle <- 1.3
k <- 3

par(mfrow = c(1, 2))
load("law82_poisson_comparison.RData")
plot(B.seq, t1.time[, k], type = "n", ylim = c(min(t1.time[, k]), max(t2.time[, k])), 
     ylab = "time in seconds", xlab = "number of bootstraps", cex.lab = cl, cex.axis = ca,
     main = "N = 82")
legend("topleft", text.col = c("blue", "orange"), 
       col = c("blue", "orange"), bty = "n", cex = cle,
       legend = c("multinomial frequencies",
                  "poisson frequencies"))
lines(B.seq, t1.time[, k], col = "blue", lwd = 2)
lines(B.seq, t2.time[, k], col = "orange")
mtext("(a)", side = 3, line = 0.5, at = 9.5e+5, cex = cp)
##
load("law_poisson_comparison.RData")
plot(B.seq, t1.time[, k], type = "n", ylim = c(min(t1.time[, k]), max(t2.time[, k])), 
     ylab = "time in seconds", xlab = "number of bootstraps", cex.lab = cl, cex.axis = ca,
     main = "N = 15")
legend("topleft", text.col = c("blue", "orange"), 
       col = c("blue", "orange"), bty = "n", cex = cle,
       legend = c("multinomial frequencies",
                  "poisson frequencies"))
lines(B.seq, t1.time[, k], col = "blue", lwd = 2)
lines(B.seq, t2.time[, k], col = "orange")
mtext("(b)", side = 3, line = 0.5, at = 9.5e+5, cex = cp)
par(mfrow = c(1, 1))


## Figure 2

cl <- 1.5
cm <- 1.5
nc <- 50

## here we are sampling 1e+5 out of the 1e+6 bootstrap replications,
## to generate the qq-plots on panels c and f (otherwise the figure
## gets to heavy.)
set.seed(123)
aux <- sample(seq(1000000), 100000, replace = FALSE) 

par(mfrow = c(2, 3))
load("law82_poisson_comparison.RData")
hist(out1, nclass = nc, xlim = c(0.5, 1), probability = TRUE, cex.main = cm, 
     main = "multinomial frequencies", xlab = "correlation", cex.lab = cl, cex.axis = ca)
mtext("(a)", side = 3, line = 0.5, at = 1, cex = cp)
hist(out2, nclass = nc, xlim = c(0.5, 1), probability = TRUE, cex.main = cm, 
     main = "Poisson frequencies", xlab = "correlation", cex.lab = cl, cex.axis = ca)
mtext("(b)", side = 3, line = 0.5, at = 1, cex = cp)
qqplot(out1[aux], out2[aux], xlab = "multinomial", ylab = "Poisson", cex.lab = cl, cex.axis = ca)
mtext("(c)", side = 3, line = 0.5, at = 0.94, cex = cp)
abline(a = 0,  b = 1)
load("law_poisson_comparison.RData")
hist(out1, nclass = nc, xlim = c(0, 1), probability = TRUE, cex.main = cm, 
     main = "multinomial frequencies", xlab = "correlation", cex.lab = cl, cex.axis = ca)
mtext("(d)", side = 3, line = 0.5, at = 1, cex = cp)
hist(out2, nclass = nc, xlim = c(0, 1), probability = TRUE, cex.main = cm, 
     main = "Poisson frequencies", xlab = "correlation", cex.lab = cl, cex.axis = ca)
mtext("(e)", side = 3, line = 0.5, at = 1, cex = cp)
qqplot(out1[aux], out2[aux], xlab = "multinomial", ylab = "Poisson", cex.lab = cl, cex.axis = ca)
mtext("(f)", side = 3, line = 0.5, at = 1, cex = cp)
abline(a = 0,  b = 1)
par(mfrow = c(1, 1))

