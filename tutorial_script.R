library(devtools)
codeUrl <- paste("https://raw.githubusercontent.com", "echaibub",
                 "VectorizedNonParametricBootstrap", "master",
                 "vectorized_bootstrap_functions.R", sep = "/")
source_url(codeUrl)


## load the American law schools data from the bootstrap R package
##
data(law) 
N <- nrow(law)
x1 <- law$LSAT
x2 <- law$GPA
x <- cbind(x1, x2)


#############################################
## Sample mean example
#############################################

B <- 1e+5
set.seed(123)
system.time( b1 <- LoopMeanBootstrap(B, x1) )
system.time( b2 <- VectorizedBootstrap(N, B, SampleBootMean, x1) )
system.time( b3 <- boot(x1, ThetaMean2, B, stype = "i")$t[, 1] )
system.time( b4 <- bootstrap(seq(N), B, ThetaMean, x1)$thetastar )


par(mfrow = c(2, 2))
hist(b1, nclass = 50, main = "for loop", xlab = "mean")
hist(b2, nclass = 50, main = "vectorized bootstrap", xlab = "mean")
hist(b3, nclass = 50, main = "R/boot", xlab = "mean")
hist(b4, nclass = 50, main = "R/bootstrap", xlab = "mean")
par(mfrow = c(1, 1))


#############################################
## Pearson's sample correlation example
#############################################

B <- 1e+5
set.seed(123)
system.time( b1 <- LoopCorBootstrap(B, x1, x2) )
system.time( b2 <- VectorizedBootstrap(N, B, SampleBootCor, x1, x2) )
system.time( b3 <- boot(x, ThetaCor2, B, stype = "i")$t[, 1] )
system.time( b4 <- bootstrap(seq(N), B, ThetaCor, x1, x2)$thetastar )

par(mfrow = c(2, 2))
hist(b1, nclass = 50, xlim = c(0, 1), main = "for loop", xlab = "correlation")
hist(b2, nclass = 50, xlim = c(0, 1), main = "vectorized bootstrap", 
     xlab = "correlation")
hist(b3, nclass = 50, xlim = c(0, 1), main = "R/boot", xlab = "correlation")
hist(b4, nclass = 50, xlim = c(0, 1), main = "R/bootstrap", xlab = "correlation")
par(mfrow = c(1, 1))


#######################################################
## Two sample problem (Welch's test statistic) example
#######################################################

## load mouse data an prepare data set for boot
##
data(mouse.t) ## treatment
data(mouse.c) ## control
groups <- c(rep(1, length(mouse.t)), rep(2, length(mouse.c)))
days <- c(mouse.t, mouse.c)
mouse <- data.frame(days, groups)

B <- 1e+5
set.seed(123)
system.time( b1 <- LoopWelchsBootstrap(B, x1 = mouse.t, x2 = mouse.c) )
system.time( b2 <- Vectorized2SampleBootstrap(length(mouse.t), 
                                              length(mouse.c), 
                                              B, SampleBootWelchs, 
                                              mouse.t, mouse.c) )
system.time( b3 <- boot(mouse, ThetaWelchs, R = B, stype = "f", 
                        strata = mouse[, 2])$t[, 1] )

par(mfrow = c(1, 3))
hist(b1, nclass = 100, probability = TRUE, xlim = c(-7, 7), 
     main = "for loop", xlab = "Welch's t statistic")
hist(b2, nclass = 100, probability = TRUE, xlim = c(-7, 7), 
     main = "vectorized", xlab = "Welch's t statistic")
hist(b3, nclass = 100, probability = TRUE, xlim = c(-7, 7), 
     main = "R/boot", xlab = "Welch's t statistic")
par(mfrow = c(1, 1))

obs <- WelchsTestStat(x1 = mouse.t, x2 = mouse.c)
obs
sum(b1 >= obs)/B
sum(b2 >= obs)/B
sum(b3 >= obs)/B
tt <- t.test(mouse.t, mouse.c, var.equal = FALSE, alternative = "greater")
tt

set.seed(123)
system.time( bb2 <- Vectorized2SampleBootstrap(length(mouse.t), length(mouse.c), 
                                               1e+6, SampleBootWelchs, mouse.t, 
                                               mouse.c) )
sum(bb2 >= obs)/1e+6

xaxis <- seq(-6, 6, length.out = 1000)
densi <- dt(xaxis, tt$parameter)
hist(bb2, nclass = 300, probability = TRUE, xlim = c(-6, 6), 
     main = "vectorized", xlab = "Welch's statistic")
lines(xaxis, densi, col = "blue", lwd = 2)
abline(v = obs, col = "red", lwd = 2)

