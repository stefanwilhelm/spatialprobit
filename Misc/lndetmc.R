# http://www.weizmann.ac.il/matlab/techdoc/ref/vander.html
# Vandermonde matrix; vander(v) = A(i,j) = v(i)^(n-j), where n = length(v)
vander <- function(v) {
 n <- length(v)
 A <- outer(1:length(v), 1:length(v), function(i, j) {v[i]^(n-j)})
 return(A)
}
#vander(v=seq(1, 3, by=0.5))
#ans =
#
#    1.0000    1.0000    1.0000    1.0000    1.0000
#    5.0625    3.3750    2.2500    1.5000    1.0000
#   16.0000    8.0000    4.0000    2.0000    1.0000
#   39.0625   15.6250    6.2500    2.5000    1.0000
#   81.0000   27.0000    9.0000    3.0000    1.0000





#% PURPOSE: computes Barry and Pace MC approximation to log det(I-rho*W)
#% -----------------------------------------------------------------------
#% USAGE: out = lndetmc(order,iter,W,rmin,rmax)
#% where:      order = # of moments u'(wsw^j)u/(u'u) to examine (default = 50)
#%              iter = how many realizations are employed (default = 30)
#%                 W = symmetric spatial weight matrix (standardized)             
#% -----------------------------------------------------------------------
#% RETURNS: out = a structure variable
#%          out.lndet = a vector of log determinants for -1 < rho < 1
#%          out.rho   = a vector of rho values associated with lndet values
#%          out.up95  = an upper 95% confidence interval on the approximation
#%          out.lo95  = a lower  95% confidence interval on the approximation
#% -----------------------------------------------------------------------
#% NOTES: only produces results for a grid of 0 < rho < 1 by default
#%        where the grid ranges by 0.01 increments
#% -----------------------------------------------------------------------
#% References: Ronald Barry and R. Kelley Pace, "A Monte Carlo Estimator
#% of the Log Determinant of Large Sparse Matrices", Linear Algebra and
#% its Applications", Volume 289, Number 1-3, 1999, pp. 41-54.
#% -----------------------------------------------------------------------
# 
# 
#% Written by Kelley Pace, 6/23/97 
#% (named fmcdetnormgen1.m in the spatial statistics toolbox )
#% Documentation modified by J. LeSage 11/99
# Ported to R by Stefan Wilhelm; renamed wsw to W
#
lndetmc <- function(order=50,iter=30,W,rmin=1e-5,rmax=1) {

n=NCOL(W)

# Exact moments from 1 to oexact
td=c(0, sum(rowSums(W%*%W))/2);
oexact=length(td);
o=order;

# Stochastic moments
mavmomi=matrix(0, o, iter);     # (o x iter), default(50 x 30)
for (j in 1:iter) {
  u=rnorm(n);        # u ~ N(0, 1), (n x 1)
  v=u;
  utu=as.numeric(t(u)%*%u);      # 1 x 1
  for (i in 1:o) {
    v=W%*%v;         # W^i*u; (n x 1)
    mavmomi[i,j]=as.numeric(n*((t(u)%*%v)/(i*utu)));
  }
}
print(mavmomi)
#TODO: mavmomi[1:oexact,]=td[,rep(1, iter)];

#averages across iterations j
avmomi=rowMeans(mavmomi);      # (o x 1)

#alpha matrix
alpha=seq(rmin,rmax, by=0.01)   # (100 x 1)
#valpha=vander(alpha);          # Vandermonde matrix; vander(v) = A(i,j) = v(i)^(n-j), where n = length(v)
valpha=outer(1:length(alpha), 1:length(alpha), function(i, j) {alpha[i]^(length(alpha)-j)})
#valphaf=fliplr(valpha);         # flip matrix left to right
valphaf=valpha[,rev(1:ncol(valpha))]  # (100 x 100)
alomat=-valphaf[,(2:(o+1))];     # (100 x 50)

#Estimated ln|I-aD| using mixture of exact, stochastic moments
#exact from 1 to oexact, stochastic from (oexact+1) to o
lndetmat=alomat %*% avmomi;      # (100 x o) x (o x 1) = 100 x 1

#standard error computations
srvs=lndetmat;
sderr=(sqrt((mean(srvs*srvs)-mean(srvs)^2)/iter));

#lower bound computation
#fbound=((n*alpha^(o+1))/((o+1)*(1-alpha+eps)));  # eps in MATLAB Maß für die interne Rechengenauigkeit
fbound=((n*alpha^(o+1))/((o+1)*(1-alpha)));

#confidendence limits, with lower biased downward (more conservative)
low95=(lndetmat-1.96*sderr-fbound);
high95=(lndetmat+1.96*sderr);

#AR parameter, lower confidence limit, estimated log-det, upper confidence limit
# confide=[alpha'  low95 lndetmat high95];

out <- list()
out$rho = alpha;
out$lndet = lndetmat;
out$up95 = high95;
out$lo95 = low95;
return(out)
}


 