library(spatialprobit)
source("sem.R")
source("rtnorm.R")

################################################################################
#
# see semp_gd.m
#
################################################################################

# W-matrix from Anselin's neigbhorhood crime data set
# n=49
anselin <- read.csv("C:/Projects/MATLAB/jplv7/data/anselin.dat", sep="",
  header=FALSE, col.names=c("crime","income","hvalue","lat","long"))
latt <- anselin[,4]
long <- anselin[,5]
W <- kNearestNeighbors(latt, long, 5); # 5 nearest neighbors weight matrix
n <- nrow(W)
sige <- 1;
k <- 3
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
rho <- 0.7

AI <- solve(I_n-rho*W)
X <- matrix(rnorm(n*k), n, k)
beta <- c(-0.5, 2, -2)
z <- X%*%beta + AI %*% rnorm(n) * sqrt(sige);
y <- as.numeric(z >= 0)

# Fit SEM model with full information z
Rprof("Anselin-sem.out")
sem_fit <- sem(z, X, W, ndraw=1000, burn.in=100, thinning=1, showProgress=TRUE)
Rprof(NULL)
summaryRprof("Anselin-sem.out")
# 17.92 Sekunden; 53.68% Zeit in "-" Methode

# Fit SEM probit model with only y=0/1 information
Rprof("Anselin-semp.out")
semp_fit <- sem_probit_mcmc(y, X, W, ndraw=1000, burn.in=100, thinning=1, showProgress=TRUE)
Rprof(NULL)
summaryRprof("Anselin-semp.out")
# 11.74 Sekunden; trotzdem 41.91% Zeit in "-" Methode

par(mfrow=c(3,2))
plot(sem_fit$B[,1],   main=expression(beta[1]), type="l", col="red",
  ylim=range(sem_fit$B[,1],semp_fit$B[,1]))     # SEM
lines(semp_fit$B[,1], type="l", col="green")    # SEM Probit

plot(sem_fit$B[,2],   main=expression(beta[2]), type="l", col="red",
  ylim=range(sem_fit$B[,2],semp_fit$B[,2]))     # SEM
lines(semp_fit$B[,2], type="l", col="green")    # SEM Probit
#lines(1:350, semp_fit2$B[1:350,2], type="l", col="blue")    # SEM Probit mit Truncated Multinormal
#lines(1:350, semp_fit3$B[1:350,2], type="l", col="orange")    # SEM Probit mit Truncated Multinormal

plot(sem_fit$B[,3],   main=expression(beta[3]), type="l", col="red", ylim=range(sem_fit$B[,3],semp_fit$B[,3]))     # SEM
lines(semp_fit$B[,3], type="l", col="green")    # SEM Probit
#lines(1:350, semp_fit2$B[1:350,3], type="l", col="blue")    # SEM Probit mit Truncated Multinormal
#lines(1:350, semp_fit3$B[1:350,3], type="l", col="orange")    # SEM Probit mit Truncated Multinormal


plot(sem_fit$B[,4],   main="sige", type="l", col="red", ylim=range(sem_fit$B[,4],semp_fit$B[,4]))     # SEM
lines(semp_fit$B[,4], type="l", col="green")   # SEM Probit

plot(sem_fit$B[,5],   main="rho", type="l", col="red", ylim=range(sem_fit$B[,5],semp_fit$B[,5]))     # SEM
lines(semp_fit$B[,5], type="l", col="green")   # SEM Probit

par(mfrow=c(3,2))

d1a <- density(sem_fit$B[,1])
d1b <- density(semp_fit$B[,1])
plot(d1a, xlim=range(d1a$x, d1b$x), col="red")
lines(d1b, col="green")
abline(v=-0.5, lty=2, col="gray40")

plot(density(sem_fit$B[,5]), col="red")
lines(density(semp_fit$B[,5]), col="green")


par(mfrow=c(3,5))
plot(semp_fit, trueparam=c(-0.5, 2, -2, sige, rho))

################################################################################
#
# Teste c_sem
#
################################################################################

tmp <- sar_lndet(ldetflag=0, W, rmin=-1, rmax=+1)
detval <- tmp$detval
detval1  <- detval[,1]  # SW: avoid multiple accesses to detval[,1]
detval2  <- detval[,2]


Rprof("c_sem.out")
for (i in 1:1000) {
  a <- c_sem(rho,z,X,beta,sige,I_n,W,detval1,detval2,rep(1,n),a1=1,a2=1)
}
Rprof(NULL)
summaryRprof("c_sem.out")

# I_n - rho * W; sehr teuer auszurechnen, weil immer über new TMat eine neue dgTMatrix
# erzeugt wird; 50% der Zeit nur in new/initialize
Rprof("S.out")
for (i in 1:1000) {
 x <- I_n - rho * W
 #x <- (-rho) * W
 #diag(x) <- diag(x) + 1
}
summaryRprof("S.out")

# Triplet representation auch nicht besser
I_n2 <- as(I_n, "dgTMatrix")
W2 <-as(W, "dgTMatrix")
Rprof("S2.out")
for (i in 1:1000) {
 x2 <- I_n2 - rho * W2
}
summaryRprof("S2.out")

# dense matrix --> rattenschnell
I_n <- as.matrix(I_n)
W <- as.matrix(W)
Rprof("S3.out")
for (i in 1:1000) {
 x <- I_n - rho * W
}
summaryRprof("S3.out")
# 0.0 Sekunden



# rechne mit (I_n - rho * W) spdep
library(spdep)
example(NY_data)
W <- as_dsTMatrix_listw(listw_NY)
for (i in 1:1000) {
  x <- as_dsCMatrix_IrW(W, rho) # funktioniert nur für symmetrische Matrizen W
}
# braucht auch ewig



