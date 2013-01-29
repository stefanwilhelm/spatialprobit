# Compute marginal effects for SAR probit
library(spatialprobit)

set.seed(0)

d <- n <- 10
I_n <- sparseMatrix(1:d, 1:d, x=1)
X <- cbind(1, seq(-2, -0.1, by=0.2), seq(3, -1.5, by=-0.5))

W <- bandSparse(d, k=-c(-1,1), diagonals=list(rep(1,d), rep(1,d)))
W <- sweep(W, 1, FUN="/", rowSums(W))

################################################################################
## compute traces for W^i
# Monte Carlo estimation of tr(W^i) for i = 0..(o-1) vs. 1..o ???
o <- 100
trW.i <- tracesWi(W, o=o, iiter=10000)
trW.i[1:5,1:5]

# Exact computation of tr(W^i) for i = 1..o
trW.i2 <- matrix(NA, n, o)
Wi <- W
trW.i2[,1] <- 1        # diag(W^0)
trW.i2[,2] <- diag(W)  # diag(W^1)
for (i in 3:o) {
  Wi <- Wi %*% W       # diag(W^(i-1)) 
  trW.i2[,i] <- diag(Wi)
}

trW.i2[1:5,1:5]
# --> Index ist um 1 verschoben 1..o

################################################################################

# parameters
beta <- c(0, 1, -1)
rho <- 0.75
beff <- c(1, -1)    # (2 x 1) without intercept

S <- I_n - rho * W
QR <- qr(S)            # class "sparseQR"
mu <- solve(QR, X %*% beta)

# compute marginal effects


rhovec <- rho^(0:(o-1)) # SW: (100 x 1)   mit [1, rho^1, rho^2 ..., rho^99], see LeSage(2009), eqn (4.145), p.115
pdfz <- dnorm(as.numeric(mu))                     # standard normal pdf phi(mu)
dd   <- sparseMatrix(i=1:n, j=1:n, x=pdfz)       # dd is diagonal matrix with pdfz as diagonal (n x n)

dir      <- as.double(t(pdfz) %*% trW.i %*% rhovec /n)  # (1 x n) * (n x o) * (o x 1)

# direct impact : dy_i / d X_ir = phi((In -  rho W)^{-1} X beta_r) * beta_r  ????
avg_direct     <- dir * beff      # (p x 1)
avg_direct
# 0.08962747 -0.08962747

# mit exact traces
dir2      <- as.double(t(pdfz) %*% trW.i2 %*% rhovec /n)  # (1 x n) * (n x o) * (o x 1)
avg_direct2     <- dir2 * beff      # (p x 1)
avg_direct2
# 0.08983675 -0.08983675


b <- as.vector(dd %*% rep(1, n))
# solve(QR, b) is (n x 1); mean(solve(QR, b)) is scalar; beff is (p x 1)
avg_total <- mean(solve(QR, b)) * beff       # (p x 1)
avg_total
# 0.2316206 -0.2316206

###
# SO ISTS RICHTIG!
# anderer Weg, total zu berechnen
Y <- solve(QR, rep(1, n))
x <- dd %*% Y * beff[1]
mean(x)
# total effects
mean(dd %*% Y) * beff
 


################################################################################

## Richtige Effekte ausrechnen
# direct impact : dy_i / d X_ir = phi((In -  rho W)^{-1} X beta_r) * beta_r
# average direct impact = average across all observations i
colMeans(apply((solve(S) %*% X[,-1]) * matrix(beff, n, 2, byrow=TRUE), 2, dnorm) * matrix(beff, n, 2, byrow=TRUE))
#[1]  0.1477278 -0.1614748
# --> Scheint falsch zu sein
r <- 1
colMeans(apply((solve(S) %*% X[,2]) * diag(beff[r], nrow=n), 2, dnorm) * matrix(beff, n, 2, byrow=TRUE))


#   direct: M_r(D) = n^{-1} tr(S_r(W))           # SW: efficient approaches available, see chapter 4, pp.114/115
#    total: M_r(T) = n^{-1} 1'_n S_r(W) 1_n      # SW: Problem: 1'_n S_r(W) 1_n ist dense!
# indirect: M_r(I) = M_r(T) - M_r(D)

# Wie sieht S_r(W) fürs Probit aus?
# LeSage (2009), equation(10.10), p.294
# S_r(W) = phi[ (I_n - rho W)^{-1} I_n xbar_r \beta_r ]  (I_n - rho W)^{-1} I_n \beta_r    (n x n)
# Diagonalelemente sind direkte Effekte = tr(S_r(W)); Average Direct Effects = 1/n tr(S_r(W))
direct <- c()
total <- c()
indirect <- c()
for (r in 1:2) {
  S_r <- apply(solve(S) %*% diag(mean(X[,-1][,r]) * beff[r], nrow=n), 2, dnorm) * (solve(S) %*% diag(beff[r], nrow=n))
  diag(S_r)
  direct[r] <- mean(diag(S_r))   # average direct effect =  M_r(D) = n^{-1} tr(S_r(W))
  # r = 1: 0.5652562
  # r = 2: -0.4909551

  # total effects: M_r(T) = n^{-1} 1'_n S_r(W) 1_n
  total[r] <- mean(S_r %*% rep(1, n))
}
direct
total
indirect <- total - direct
indirect

# --> direct effects weichen stark ab zwischen
# (a) der Näherung über tr(W^i)
# (b) equation(10.10)