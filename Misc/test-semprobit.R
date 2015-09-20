library(spatialprobit)
source("G:/R/spatialprobit/R/rtnorm.R")

# SEM Probit zu Fuﬂ: UNtersuche nur mal das Sampling z_i | z_{-i}, beta, rho, sige

### Model Setup
n <- 5
X <- matrix(cbind(1, rnorm(n)), ncol=2)
rho <- 0.3
beta <- c(1, -1)
sige <- 1.2
I_n  <- sparseMatrix(i=1:n, j=1:n, x=1)
W <- bandSparse(5, k=c(-1, 1), diag=list(rep(1, 5), rep(1, 5)))
W <- W / rowSums(W)
ones <- rep(1, n)
W2diag <- diag(t(W)%*%W)

S    <- (I_n-rho*W)
AI   <- solve(I_n-rho*W)
z    <- X%*%beta + AI %*% rnorm(n) * sqrt(sige)
y    <- as.numeric(z >= 0)

ind0 <- which(y == 0)      # 0 obs.
ind1 <- which(y == 1)      # 1 obs.

# truncation points for z, depend only on y, can be precalculated
lower <- ifelse(y > 0, 0,  -Inf)
upper <- ifelse(y > 0, Inf,   0)



### Estimate z_i | z_{-i}, beta, rho, sige

H <- 1/sige * t(S) %*% S
diag(H)
1/sige * (diag(I_n) + rho^2*W2diag)

# conditional variance  z_i | z_{-i}
dsig <- 1/sige * (ones + rho * rho * W2diag)  # SW: Sollte es nicht ones + rho * rho * W2diag sein?
zvar <- ones/dsig;            # conditional variances for each z_i | z = H_{ii}^{-1} = 1/diag(H) (n x 1)
                              # TODO: Was passiert im Fall dsig < 0 --> zvar < 0 (negative variance) !!!
all.equal(zvar, 1 / diag(H))    # TRUE zvar = diag(H)^{-1}

# conditional mean  z_i | z_{-i}
mu <- X %*% beta
zmu <- z - mu
A  <- (1/sige)* S %*% zmu     # a vector (n x 1)
B2  <- t(S) %*% A             # B2 <- (1/sige) * t(S) %*% S %*% zmu
Cz <- zmu - zvar*B2           # Cz = (z - mu) - diag(H)^{-1} * H * (z - mu)
zm <- mu + Cz;                # mu + (z-mu) - zvar * t(S)[ (1/sige) * S * (z - mu) ]  =  mu + (z-mu)
# conditional mean based on H
zm2 <- mu - 1 / diag(H) * H %*% zmu  + zmu
all.equal(zm, zm2)          # TRUE


# sampled based on univariate truncated normal
z1 <- matrix(NA, 100, n)
# sampled based on multivariate truncated normal
z2 <- matrix(NA, 100, n)
for (i in 1:100) {
  z1[i,ind0] <- rtnorm(mu=zm[ind0], sd=sqrt(zvar[ind0]), a=-Inf, b=0)
  z1[i,ind1] <- rtnorm(mu=zm[ind1], sd=sqrt(zvar[ind1]), a=0, b=Inf)
  
  
  z2[i,] <- rtmvnorm.sparseMatrix(n=1, mean=mu, H=H, lower=lower, upper=upper,
                      start.value=as.double(z), burn.in=10)
}

plot(ecdf(z1[,1]))
lines(ecdf(z2[,1]), col="red")

plot(ecdf(z1[,2]))
lines(ecdf(z2[,2]), col="red")




