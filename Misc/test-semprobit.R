library(Matrix)
library(tmvtnorm)


# Teste SEM model: (n=2)
# z = X beta + (I_n-rho*W)*eps, eps ~ N(0, sige^2 * I_n)
n <- 2
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
rho <- 0.5
W <- matrix(c(0.5, 0.5,
              0.5, 0.5), 2, 2)
sige <- 2
A <- (I_n - rho * W)
AI <- solve(A)
X <- matrix(1, 2, 1)
beta <- 0
z <- X%*%beta + solve(A) %*% (sige*rnorm(n=2, sd=1))
#--> Modell führt zu einer Korrelation der Residuen von 0.6!

# simuliere N Wdh. aus dem Modell z = Xbeta + (I_n-rho*W)^(-1)*eps, eps ~ N(0, sige^2 * I_n)
N <- 10000
eps <- matrix(rnorm(n*N), N, n)   # (N x n)
Z <- as.matrix(sige * eps %*% AI)
plot(Z[,1], Z[,2])
cor(Z)            # cor=0.5748

# simuliere direkt aus der Verteilung z ~ N(Xbeta, sige^2*AI'AI)
Sigma <- as.matrix(sige^2*t(AI)%*%AI)
Z2 <- rmvnorm(n=10000, mean=X%*%beta, sigma=Sigma)
points(Z2, col="red")
cor(Z2)

# Teste Truncation:  Y = Z > 0
Y <- as.matrix(Z>0)
Y2 <- as.matrix(Z2>0)

# Teste Simulation von z | y, beta, rho, sige ~ TMVN(Xbeta, sige^2*AI'AI, lower, upper)
# Bsp. z[1] < 0, z[2] > 0 --> y[1]=0, y[2]=1
lower<-c(-Inf,0)
upper<-c(0, Inf)
Z3 <- rtmvnorm(n=10000, mean=as.numeric(X%*%beta), sigma=Sigma, lower=lower, upper=upper)

par(mfrow=c(2,2))
plot(density(Z3[,1]))
ind <- Y2[,1]==0&Y2[,2]==1
lines(density(Z2[ind,1]), col="red")
# --> passt