# test SAR Tobit model using Koop example; see sart_gd.m from LeSage
library(spatialprobit)
n <- 300

# generate explanatory variable
a <- 0
b <- 1
x <- a+(b-a)*rnorm(n)

X <- cbind(rep(1, n), x)

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

y <- S %*% (X%*%beta) + S %*% e
ysave <- y

ind <- which(y <=0)
y[ind] <- 0
nobsc <- length(ind)

# Fit SAR (with complete information) as double check if beta, sige are estimated correctly
Rprof("sar2.out")
fit_sar <- sartobit(ysave,X,W,ndraw=200,burn.in=100, showProgress=FALSE)
Rprof(NULL)
summaryRprof("sar2.out")
par(mfrow=c(2, 4))
plot(fit_sar, which=c(1, 2), trueparam=params)
# --> Passt fürs SAR. Kann man so durchgehen lassen
# ca. 8.68 ohne ProgressBar; 46.08% davon in "-" für CSparse

# SAR Tobit
Rprof("sartobit.out")
fit_sart <- sartobit(y,X,W,ndraw=200,burn.in=100, showProgress=TRUE)
Rprof(NULL)
summaryRprof("sartobit.out")
summary(fit_sart)

par(mfrow=c(3, 4))
plot(fit_sart, which=c(1, 2), trueparam=params)
