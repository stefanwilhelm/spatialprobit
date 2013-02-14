library(Matrix)
library(spdep)

################################################################################
#
# QR decomposition for dense Matrix
#
################################################################################

I_n <- diag(3)
W <- matrix(c(0, 1, 0,
              1, 0, 1,
              0, 1, 0), 3, 3)
rho <- 0.75
S <- (I_n - rho * W)
qrS <- qr(S)     # S = QR

R <- qr.R(qrS)   # upper triangular matrix
Q <- qr.Q(qrS)   # orthogonal matrix Q'Q = I

Q %*% R  # --> should give S

# S^(-1) = R^(-1) * Q^T
system.time(for (i in 1:10000) solve(S))
system.time(for (i in 1:10000) solve(R) %*% t(Q))
# 0.78 vs. 1.10
all.equal(solve(S), solve(R) %*% t(Q))

# S^(-1) %*% 1_n
solve(S) %*% 1

# Diagonalmatrix mit phi(z)
D <- diag(c(0.2, 0.3, 0.5))

# Koeffizient beta
b <- 2

# gesucht (a)
colSums(D %*% solve(S) * b)

# (b) gleichwertig ohne Invertierung von S,
# nur über die QR-Zerlegung von S die wir ohnehin schon haben
D %*% rep(1, n) * b
solve(qr(S), D %*% rep(1, n) * b)
solve(qr(S), D %*% rep(1, n)) * b


################################################################################
#
# QR decomposition for sparse Matrix
#
################################################################################

n <- 3
I_n <- sparseMatrix(1:n, 1:n, x=1)
W <- Matrix(c(0, 1, 0,
              1, 0, 1,
              0, 1, 0), 3, 3, sparse=TRUE)
              
S <- (I_n - rho * W)

qrS <- qr(S)

qrS@R
