library(spatialprobit)

################################################################################
#
## optimize tr(W^i)
#
################################################################################

# Vermeide doppelte Schleifen für die Berechnung von u'(W^i)u...

# estimated tr(W^i) for i=1,...,100
# see LeSage (2009), chapter 4, p.97/98
#
# "The expectation of the quadratic form u'Au equals tr(A).
#  Since u_i^2 follows a chi^2 distribution with one degree of freedom."
# Pre-calculate traces tr(W^i) i=1,...,100 for the x-impacts calculations
#
# @param W spatial weight matrix (n x n)
# @param o highest order of W^i = W^o
# @param iiter number of MCMC iterations (we run this 50 times)
# @return (n x o) matrix with tr(W^i) in each column, for i=1..,o
tracesWi <- function(W, o=100, iiter=50) {
  n <- nrow(W)
  trW_i <- matrix( data=0, nrow=n, ncol=o )   # n x 100
  for( iii in 1:iiter ){
    u        <- rnorm(n)         #u    -> n x 1      # randn() aus Paket pracma
    wumat    <- matrix(0, n, o)  #wumat-> n x 100
    wumat[,1]<-u             
    wu       <-u                 #wu    -> n x 1
    for(i in 2:o ){
      wu         <- W %*% wu     # (n x n) * (n x 1) = (n x 1) matrix
      wumat[ ,i] <- as.double(wu)
    }
    trW_i        <- trW_i + u * wumat  # n x 100    u' W^i u  (i = 1..100)
  }
  trW_i <- trW_i / iiter
}

tracesWi_2 <- function(W, o=100, iiter=50) {
  n <- nrow(W)
  trW_i <- matrix( data=0, nrow=n, ncol=o )   # n x o
  u  <- matrix(rnorm(n * iiter), nrow = n, ncol = iiter)   # (n x iiter)
  xx <- u
  trW_i[,1] <- apply(u * as.matrix(xx), 1, sum)
  for(i in 2:o ){
    xx <- W %*% xx  # (n x iter)
    trW_i[,i] <- apply(u * as.matrix(xx), 1, sum)  # u'(W^i)u; sum across all iterations
  }
  trW_i <- trW_i / iiter
  return(trW_i)
}

d <- n <- 10
W <- bandSparse(d, k=-c(-1,1), diagonals=list(rep(1,d), rep(1,d)))
W <- sweep(W, 1, FUN="/", rowSums(W))

set.seed(1)
tr1 <- tracesWi(W, o=10, iiter=20)
set.seed(1)
tr2 <- tracesWi_2(W, o=10, iiter=20)
all.equal(tr1, tr2)

Rprof("tracesWi.out")
set.seed(1)
for (i in 1:10)
tr1 <- tracesWi(W, o=100, iiter=50)
Rprof(NULL)
summaryRprof("tracesWi.out")

Rprof("tracesWi2.out")
set.seed(1)
for (i in 1:10)
tr2 <- tracesWi_2(W, o=100, iiter=50)
Rprof(NULL)
summaryRprof("tracesWi2.out")
# --> Optimierte Methode tracesWi_2 braucht nur noch 1/10 der Zeit von tracesWi!
# --> vgl. mcdet_setup
# trWi: 10 Sekunden
# trWi2: 0.5 Sekunden
all.equal(tr1, tr2)
