# Example from LeSage/Pace (2009), section 10.3.1, p. 302-304
#
# n = 1000 obs
# x_i ~ U(a, 1)
# eps_i ~ N(0, 0.5), sige=0.5
# beta = c(0, 2)
# rho = 0.7
# y = (I_n - rho * W)^{-1}(X beta + eps)
# 6 nearest neighbors
# Value of a is not stated in book! Assuming a=-1 which gives approx. 50% censoring
# ndraw = 1000, 200 burn-in

library(spatialprobit)

a <- -1   # control degree of censored observation
n <- 1000
rho <- 0.7
beta <- c(0, 2)
sige <- 0.5
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
x <- runif(n, a, 1)
X <- cbind(1, x)
eps <- rnorm(n, sd=sqrt(sige))
param <- c(beta, sige, rho)

# random locational coordinates and 6 nearest neighbors
lat <- rnorm(n)
long <- rnorm(n)
W <- kNearestNeighbors(lat, long, k=6)

y <- as.double(solve(I_n - rho * W) %*% (X %*% beta + eps))
table(y > 0)

# full information
ysave <- y

# set negative values to zero to reflect sample truncation
ind <- which(y <=0)
y[ind] <- 0

# Fit SAR (with complete information) as double check if beta, sige and rho are estimated correctly
Rprof("sar.out")
fit_sar <- sartobit(ysave ~ X-1, W,ndraw=1000,burn.in=200, showProgress=FALSE)
Rprof(NULL)
summaryRprof("sar.out")
# 3 Sekunden
summary(fit_sar)

#----MCMC spatial autoregressive Tobit model ----
#Execution time  = 3.379 secs
#
#N draws         =   1000, N omit (burn-in)=    200
#N observations  =   1000, K covariates    =      2
## censored values =      0, # observed values =   1000
#Min rho         = -1.000, Max rho         =  1.000
#--------------------------------------------------
#
#     Estimate Std. Dev p-level t-value Pr(>|z|)
#      0.01818  0.02276 0.20600   0.799    0.425
#x     2.02408  0.03893 0.00000  51.988   <2e-16 ***
#sige  0.54311  0.02501 0.00000  21.714   <2e-16 ***
#rho   0.69440  0.01767 0.00000  39.300   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

par(mfrow=c(3, 4))
plot(fit_sar, which=c(1, 2, 3), trueparam=param)

# SAR Tobit
Rprof("sartobit.out")
fit_sartobit <- sartobit(y ~ x,W,ndraw=1000,burn.in=200, showProgress=TRUE)
Rprof(NULL)
summaryRprof("sartobit.out")
summary(fit_sartobit)
# 19 Sekunden mit Progressbar; 15.574 Sekunden ohne

#----MCMC spatial autoregressive Tobit model ----
#Execution time  = 16.549 secs
#
#N draws         =   1000, N omit (burn-in)=    200
#N observations  =   1000, K covariates    =      2
## censored values =    471, # observed values =    529
#Min rho         = -1.000, Max rho         =  1.000
#--------------------------------------------------
#
#     Estimate Std. Dev p-level t-value Pr(>|z|)
#      0.02194  0.02715 0.21900   0.808    0.419
#x     2.06971  0.06969 0.00000  29.700   <2e-16 ***
#sige  0.57208  0.03919 0.00000  14.597   <2e-16 ***
#rho   0.67536  0.02018 0.00000  33.475   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

par(mfrow=c(3, 4))
plot(fit_sartobit, which=c(1, 2, 3), trueparam=param)


par(mfrow=c(2,2))
ylim1 <- range(fit_sar$B[,1], fit_sartobit$B[,1])
plot(fit_sar$B[,1], type="l", ylim=ylim1, col="red")
lines(fit_sartobit$B[,1], col="green")

ylim2 <- range(fit_sar$B[,2], fit_sartobit$B[,2])
plot(fit_sar$B[,2], type="l", ylim=ylim2, col="red")
lines(fit_sartobit$B[,2], col="green")

ylim3 <- range(fit_sar$B[,3], fit_sartobit$B[,3])
plot(fit_sar$B[,3], type="l", ylim=ylim3, col="red")
lines(fit_sartobit$B[,3], col="green")





