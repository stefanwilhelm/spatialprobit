require(pracma)

# PURPOSE: constructs a row-stochastic nearest neighbor spatial weight matrix
#          asymmetric, but row-sums are unity, based on m neighbors
# --------------------------------------------------------
# USAGE: W = make_neighborsw(xc,yc,nn)
#       where: 
#             xc = x-coordinate for each obs (nobs x 1)
#             yc = y-coordinate for each obs (nobs x 1)
#             nn = # of nearest neighbors to be used
# --------------------------------------------------------
# RETURNS: W an (nobs x nobs) spatial weight matrix based on nn
#          nearest neighbors (a sparse matrix)
# --------------------------------------------------------
# NOTES: 
#        W takes a form such that: W*y would produce a vector
#        consisting of the values of y for the nn nearest neighbors
#        for each observation i in the (nobs x 1) vector y
# To construct a weight matrix based on 4 nearest neighbors
# W4 = make_neighborsw(xc,yc,4);
#   ---> This function will is similar to make_nnw, but uses less
#        memory and takes more time. If you run out of memory using
#        make_nnw, try this function
# --------------------------------------------------------
# SEE ALSO: find_neighbors(), find_nn(), make_nnw()
# --------------------------------------------------------
# written by:
# James P. LeSage, 5/2002
# updated 1/2003
# Dept of Economics
# University of Toledo
# 2801 W. Bancroft St,
# Toledo, OH 43606
# jlesage@spatial-econometrics.com
#
#Adapted by:
#Miguel Godinho de Matos
#Dept of Engineering & Public Policy
#Carnegie Mellon University
make_neighborsw<-function(xc,yc,m){
  n      <- nrow(xc)
  nnlist <- find_neighbors(xc,yc,m)
  #convert the list into a row-standardized spatial weight matrix
  rowseqs <- (1:n) 
  vals1   <- ones(n,1)*(1/m)
  vals0   <- zeros(n,1);
  for( i in 1:m ){
    colseqs     <- nnlist[ , i]
    ind_to_keep <- colseqs > 0 
    z1          <- cbind( rowseqs, colseqs, vals1 )
    z1          <- z1[ ind_to_keep, ]
    #z2          <- cbind( rowseqs, rowseqs, vals0 )
    z           <- rbind( z1, c(n,n,0) ) #z2 sets the size, but z1 might be wrong
    if( i == 1 ){
      W <- spconvert( z )
    }else{
      W <- W + spconvert( z )
    }
  }
  return( W )
}

# PURPOSE: finds observations containing m nearest euclidean distance-based neighbors,
#          (slow but low memory version) returns an index to these neighboring observations
# --------------------------------------------------------
# USAGE: nnindex = find_neighbors(xc,yc,m)
#       where: 
#             xc = x-coordinate for each obs (nobs x 1)
#             yc = y-coordinate for each obs (nobs x 1)
#             m  = # of nearest neighbors to be found
# --------------------------------------------------------
# RETURNS: an (nobs x m) matrix of indices to the m neighbors
# --------------------------------------------------------
# NOTES: nnindex takes a form such that: ind = nnindex(i,:)';
#        y(ind,1) would pull out the m nearest neighbor observations to
#        y(i,1), and y(ind,1)/m would represent an avg of these
#   ---> This function will is similar to find_nn, but uses less
#        memory and takes more time. If you run out of memory using
#        find_nn, try this function
# --------------------------------------------------------
# SEE ALSO: find_nn, make_neighborsw, make_nnw, make_xyw
# --------------------------------------------------------

# written by:
# James P. LeSage, 12/2001
# modified 1/2003
# Dept of Economics
# University of Toledo
# 2801 W. Bancroft St,
# Toledo, OH 43606
# jlesage@spatial-econometrics.com
# NOTE: this is a fast approach, but requires a lot of RAM memory

#Adapted by:
#Miguel Godinho de Matos
#Dept of Engineering & Public Policy
#Carnegie Mellon University

find_neighbors <- function(xc,yc,m){
  dimx <- dim( xc )
  dimy <- dim( yc )
  n    <- dimx[ 1 ]
  
  if( dimx[1] != dimy[1] ){
    print('find_neighbors: xc,yc inputs must be same size')
    return(NULL)
  }else if( dimx[2] != 1 || dimy[2] != 1){
    print('find_neighbors: expectin n x 1 arrays for xc and xy')
    return(NULL)
  }
  
  nnlist <- zeros(n,m)
  for( i in 1:n ){
    xi   <- xc[i,1]
    yi   <- yc[i,1]
    dist <- (xc - xi * ones(n,1))^2 + (yc - yi*ones(n,1))^2
    tmp  <- sort( dist , index.return=TRUE)
    xds  <- tmp$x
    xind <- tmp$ix
    nnlist[i,] <- xind[2:(m+1)] #the closest neighbor will be itself ...
  }
  return( nnlist )
}
