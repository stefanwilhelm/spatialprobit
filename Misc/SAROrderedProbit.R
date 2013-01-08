# SAR Ordered Probit / Ordered spatial probit model
#
# see LeSage (2009), section 10.2
#
# model:
# (1)  z = rho * W  * z + X beta + eps
# (2) y_i can take J alternatives for
#     y_i = j, if phi_{j-1} <= z <= phi_j
#
# SAR probit is special case with J=2 and phi=c(-Inf, 0, Inf) is J+1 vector with phi_0 = -Inf and phi_J = +Inf
# Model paramaters to be estimated:
# beta, rho and vector phi (J-2 values)

library(tmvtnorm)
library(spatialprobit)

################################################################################
#
# Example with J = 3
#
################################################################################

# set up a model like in SAR probit
J <- 4   # ordered alternatives j=1, 2, 3, 4 --> 2 cutoff-points to be estimated phi_2, phi_3
phi <- c(-Inf, 0,  +1, +2, Inf)    # phi_0,...,phi_j, vector of length (J+1)
# phi_1 = 0 is a identification restriction

# generate random samples from true model
n <- 400             # number of items
k <- 3               # 3 beta parameters
beta <- c(0, -1, 1)  # true model parameters k=3 beta=(beta1,beta2,beta3)
rho <- 0.75
# design matrix with two standard normal variates as "coordinates"
X <- cbind(intercept=1, x=rnorm(n), y=rnorm(n))

# identity matrix I_n
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)

# build spatial weight matrix W from coordinates in X
W <- kNearestNeighbors(x=rnorm(n), y=rnorm(n), k=6)

# create samples from epsilon using independence of distributions (rnorm()) to avoid dense matrix I_n
eps <- rnorm(n=n, mean=0, sd=1)
z   <- solve(qr(I_n - rho * W), X %*% beta + eps)

# ordered variable y:
# y_i = 1 for phi_0 < z <= phi_1; -Inf < z <= 0
# y_i = 2 for phi_1 < z <= phi_2
# y_i = 3 for phi_2 < z <= phi_3
# y_i = 3 for phi_3 < z <= phi_4

# y in {1, 2, 3}, 
y   <- cut(as.double(z), breaks=phi, labels=FALSE, ordered_result = TRUE)
table(y)

################################################################################
S <- I_n - rho * W
H <- t(S) %*% S

# immutable matrices
c    <- matrix(0, k, 1)
Tinv <- matrix(0, k, k)       # uninformative prior

tX <- t(X)                    # X'               # k x n
xpx  <- t(X) %*% X            # (X'X)            # k x k
xpxI <- solve(xpx)            # (X'X)^{-1}       # k x k
xxpxI <- X %*% xpxI           # X(X'X)^(-1)     # n x k (better, compromise)
AA    <- solve(xpx + Tinv)    # (X'X + T^{-1})^{-1}

# MCMC
ndraw <- 1000
params <- k + 1 + (J+1)            # 7 parameters rho, beta (k), phi (J + 1), aber nur J-2 zu schätzen
# matrix to store the beta + rho parameters for each iteration/draw
B <- matrix(NA, ndraw, params)
colnames(B) <- c("rho", paste("beta_", 1:k, sep=""), paste("phi_", 0:J, sep=""))


phi <- c(-Inf, 0, 0.5, 2.5, Inf)   # start vector

# MCMC loop  
for (i in 1:ndraw) {

  # 1. sample from z | rho, beta, y using precision matrix H
  mu <- solve(qr(S), X %*% beta)
  
  # determine lower and upper bounds for z depending on the value of y and phi
  # eqn (10.14), p.299
  
  #lower <- rep(NA, n)
  #upper <- rep(NA, n)
  #for (j in 1:J) {
  #  ind <- y == j
  #  lower[ind] <- phi[j]
  #  upper[ind] <- phi[j+1]
  #}
  
  # besser: y = 1..J
  lower <- phi[y]
  upper <- phi[y+1]
  
  # see cbind(y, lower, upper)

  # see LeSage (2009) for choice of burn-in size, often m=5 or m=10 is used!
  # we can also use m=1 together with start.value=z, see LeSage (2009), section 10.1.5
  z <- as.double(rtmvnorm.sparseMatrix(n=1, mean=mu, H=H,
      lower=lower, upper=upper, burn.in=10))

  # 2. sample from beta | rho, z, y
  c <- AA  %*% (tX %*% S %*% z + Tinv %*% c)
  T <- AA   # no update basically on T, TODO: check this
  beta <- as.vector(rmvnorm(n=1, mean=c, sigma=T))

  # 3. sample from rho | beta, z
  # keep rho fixed for now
  rho <- 0.75
  
  # 4. determine bounds/cut-points p(phi_j | phi_{-j}, z, y, beta) for j = 2,...,J-1
  for (j in 2:(J-1)) {
    phi.lower <- max(max(z[y == j]),     phi[j-1+1])   # \bar{phi}_{j-1}, SW: +1 is needed as our vector index starts with 1
    phi.upper <- min(min(z[y == j + 1]), phi[j+1+1])   # \bar{phi}_{j+1}
    
    # Sample phi_{j | phi_{-j}, z, y, beta)
    phi[j + 1]   <- runif(n=1, min=phi.lower, max=phi.upper)
  }

  # update S and H
  S <- I_n - rho * W
  H <- t(S) %*% S      # H = S'S  / SW: crossprod(S) does not seem to work!
  
  # save parameters in this MCMC round
  B[i, ] <- c(beta, rho, phi)
}

par(mfrow=c(2,2))
plot(B[,7], type="l")        # path of phi_2
abline(h=1, lty=2, col="red")

plot(B[,8], type="l")        # path of phi_3
abline(h=2, lty=2, col="red")



