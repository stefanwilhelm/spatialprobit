library(spatialprobit)
library(tmvtnorm)
library(Matrix)

n <- d <- 300
m <- 3
W <- sparseMatrix(i=rep(1:d, each=m),
  j=replicate(d, sample(x=1:d, size=m, replace=FALSE)), x=1/m, dims=c(d, d))
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
sige <- 1.2                                # Varianz!
eps <- rnorm(n=n, mean=0, sd=sqrt(sige))   # Normierung!!!
rho <- 0.75
X   <- cbind(x1=1, x2=runif(n=n, -2, 2))
beta <- c(-0.2, 0.5)

z <- as.vector(X %*% beta + solve(I_n - rho * W) %*% eps)     # SEM model
y <- as.numeric(z >= 0)

# Faktor c
c <- 2
z2 <- as.vector(c * X %*% beta + c * solve(I_n - rho * W) %*% eps)     # SEM model
y2 <- as.numeric(z2 >= 0)
all.equal(y, y2)

mu <- X %*% beta

# truncation points for z, depend only on y, can be precalculated
lower <- ifelse(y > 0, 0,  -Inf)
upper <- ifelse(y > 0, Inf,   0)

S <- I_n - rho * W
H <- t(S) %*% S            # precision matrix H for beta | rho, z, y

znew <- as.double(rtmvnorm.sparseMatrix(n=1, mean=mu, H=H, lower=lower, upper=upper, burn.in=10))

# stürzt in rtmvnorm.sparseMatrix ab, wenn I_n/W dense "matrix" sind.
# dann ist auch H eine dense matrix. TODO: Fehler abfangen!
Rprof("SEM-probit.out")
fit <- sem_probit_mcmc(y, X, W, ndraw=200, burn.in=100, thinning=1, showProgress=TRUE)
# 1. ndraw=200, burn.in=100 : 14 Sekunden
Rprof(NULL)
summaryRprof("SEM-probit.out")

# now with dense matrix W,I_n
Rprof("SEM-probit-dense.out")
fit2 <- sem_probit_mcmc(y, X, W=as.matrix(W), ndraw=200, burn.in=100, thinning=1, showProgress=TRUE)
# 1. ndraw=200, burn.in=100 : 14 Sekunden
Rprof(NULL)
summaryRprof("SEM-probit-dense.out")


summary(fit)
par(mfrow=c(4,3))
plot(fit)
head(fit$B, 10)

tmp <- sar_lndet(ldetflag=0, W, rmin=-1, rmax=+1)
detval <- tmp$detval
detval1 <- detval[,1]
detval2 <- detval[,2]
ones <- rep(1,n)
Rprof("c_sem.out")
for (i in 1:1000) {
  a <- c_sem(rho,z,X,beta,sige,I_n,W,detval1,detval2,ones,a1=1,a2=1)
}
Rprof(NULL)
summaryRprof("c_sem.out")

I2 <- as(I_n,"TsparseMatrix")
W2 <- as(W,"TsparseMatrix")
system.time(for (i in 1:1000) A <- I_n - 0.75 * W)   # 3.62 Warum dauert das so lange?
system.time(for (i in 1:1000) A <- I2 - 0.75 * W2)   # 6.05
I3 <- as.matrix(I_n)
W3 <- as.matrix(W)
system.time(for (i in 1:1000) A <- I3 - 0.75 * W3)   # 2.40