# test SAR Tobit model using Koop example; see sart_gd.m from LeSage
library(spatialprobit)
n <- 300

# generate explanatory variable
a <- -1
b <- 1
x <- a+(b-a)*rnorm(n)

X <- cbind(intercept=rep(1, n), x2=x)

beta <- c(-1, 2)
sigma1 <- 2
rho <- 0.7
params <- c(beta, sigma1, rho)

latt <- rnorm(n)
long <- rnorm(n)
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
W <- kNearestNeighbors(latt,long,6)


S <- solve(I_n - rho*W)
e <- rnorm(n)*sqrt(sigma1)

y <- as.numeric(S %*% (X%*%beta) + S %*% e)
ysave <- as.numeric(y)

ind <- which(y <=0)
y[ind] <- 0
nobsc <- length(ind)   # number of censored observations

# Fit SAR (with complete information) as double check if beta, sige are estimated correctly
Rprof("sar.out")
fit1 <- sartobit(ysave ~ x, W,ndraw=1000, burn.in=200, showProgress=FALSE)
Rprof(NULL)
summaryRprof("sar.out")

par(mfrow=c(2, 4))
plot(fit1, which=c(1, 2), trueparam=params)
# --> Passt fürs SAR. Kann man so durchgehen lassen
# ca. 8.68 ohne ProgressBar; 46.08% davon in "-" für CSparse

# SAR Tobit with approx. 50% censored observations
#Rprof("sartobit.out")
fit2 <- sartobit(y ~ x, W,ndraw=1000,burn.in=200, showProgress=TRUE)
#Rprof(NULL)
#summaryRprof("sartobit.out")
summary(fit2)

par(mfrow=c(2,2))
plot(fit1$B[,1], type="l", ylim=range(fit1$B[,1], fit2$B[,1]), col="red")
lines(fit2$B[,1], col="green")

plot(fit1$B[,2], type="l", col="red")
lines(fit2$B[,2], col="green")



par(mfrow=c(3, 4))
plot(fit_sart, which=c(1, 2), trueparam=params)
