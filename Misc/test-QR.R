# How to solve Sx=1 using the QR-decomposition of S?
n <- 200
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
rho <- 0.75
W <- kNearestNeighbors(x=rnorm(n), y=rnorm(n))

S <- (I_n - rho * W)
QR <- qr(S)

Rprof("way1.out")
for (i in 1:500)
y <- solve(QR, rep(1, n))
Rprof(NULL)
summaryRprof("way1.out")
# n=100; i=500 --> 1.16 Sek
# n=200; i=500 --> 1.2 Sek

Rprof("way2.out")
for (i in 1:500)
y2 <- qr.coef(QR, rep(1, n))
Rprof(NULL)
summaryRprof("way2.out")
# n=100; i=500 --> 0.4 Sekunden
# n=200; i=500 --> 0.4

all(y == y2)
# TRUE

# Pre-select specific method does not seem to help much!
selectMethod("solve", signature=c("sparseQR", "numeric"))
m <- selectMethod("qr.coef", signature=c("sparseQR", "numeric"))
Rprof("way3.out")
for (i in 1:500)
y2 <- m(QR, rep(1, n))
Rprof(NULL)
summaryRprof("way3.out")
