# problem: using a generic method like qr() requires a lot of time in "searching"
# the specific method. If we find the specific method once, is it faster then?

library(spdep)

f <- findMethod("qr", "dgCMatrix")
f <- getMethod("qr", "dgCMatrix")
str(f)

# SE_Classic_Setup --> SE_setup_intern --> LU_setup
n <- 100
W <- bandSparse(n=n, k=c(-1,1), diagonals=list(c(1, rep(0.5, n-2)), c(rep(0.5, n-2), 1)))
I <- as_dsCMatrix_I(100)
# --> do_ldet rechnet nur 100 Gridpoints von ln|I - rho W| aus und interpoliert dann auf 1000/2000 Gridpoints
# dauert ca. 1.2 Sekunden für 100 Gridpoints, also auch ziemlich langsam
# --> Wir müssen ca. 1000x die Matrix (I-rho W) ausrechnen.

rho <- 0.7
I <- sparseMatrix(i=1:n,j=1:n,x=1)
S <- (I - rho * W)

# use the pre-selected method for qr("dgCMatrix")
Rprof("qr1.out")
system.time(
 for (i in 1:10000) QR <- f(S)
)
Rprof(NULL)
summaryRprof("qr1.out")

# n=10; 1000x = 0.22 Sekunden
# n=100; 1000x = 0.28 Sekunden

# use the generic method qr()
Rprof("qr2.out")
system.time(
 for (i in 1:10000) QR <- qr(S)
)
Rprof(NULL)
summaryRprof("qr2.out")
#    user  system elapsed
#   0.26    0.00    0.27
# n=100; 1000x = 0.33 Sekunden
# --> etwas lansamer
# --> erklärt beides nicht wirklich, warum soviel Zeit in "standardGeneric" verbracht wird

