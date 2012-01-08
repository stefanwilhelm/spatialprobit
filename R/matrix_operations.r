require(Matrix)
require(pracma)

# PURPOSE: Creates an nxn sparse identity matrix 
#---------------------------------------------------
# USAGE: result = speye(n)
#---------------------------------------------------
# RETURNS: a sparse n x n identity matrix
# --------------------------------------------------
speye <- function(n){
  In <- Matrix(data=0,ncol=n,nrow=n,sparse=T)
  diag(In) <- 1
  return(In)
}

# PURPOSE: converts a three column matrix in a 
#          sparse matrix
#---------------------------------------------------
# USAGE: result = spconvert(X)
# where: X is a three colum matrix that should be 
#        interpreted as [ nrow ncol value ]
#---------------------------------------------------
# RETURNS: A  martix where W[nrow,ncol] = value
#          as defined by each row of the X matrix
# --------------------------------------------------
spconvert <- function(X){
  n      <- max(X[,1])
  m      <- max(X[,2])
  Xrlist <- which( X[,3] != 0 )
  result <- Matrix(data=0,nrow=n,ncol=m,sparse=T)
  for(Xr in Xrlist ){
    result[ X[Xr,1], X[Xr,2] ] <- X[Xr, 3]
  }
  return(result)
}

# PURPOSE: performs matrix division even if matrices
#          are not of the same dimension, but are row or
#          column compatible
#---------------------------------------------------
# USAGE: result = matdiv(x,y)
# where:    x,y = two matrices (not of the same dimension
#                 but are row or column compatible)
#---------------------------------------------------
# RETURNS: result = x / y where x and y are row or column compatible
# --------------------------------------------------
# written by:
# James P. LeSage, Dept of Economics
# Unioutersity of Toledo
# 2801 W. Bancroft St,
# Toledo, OH 43606
# jpl@jpl.econ.utoledo.edu
matdiv <- function(x,y){
  rx <- nrow(x); cx <- ncol(x)
  ry <- nrow(y); cy <- ncol(y)
  if( (cx == cy) && (rx == ry) ){
    results <- x / y
  }else if( (cx == cy) && (rx == 1) ){
    results <- y / repmat(x,ry,1) 
  }else if( (cx == cy) && (ry == 1) ){
    results <- x / repmat(y,rx,1)
  }else if( (rx == ry) && (cx == 1) ){
    results <- y / repmat(x,1,cy)
  }else if( (rx == ry) && (cy == 1) ){
    results <- x / repmat(y,1,cx) 
  }else{
    print('matdiv: non-conformable in row or column dimension')
    results <- NULL
  }
  return(results)
}

# PURPOSE: performs matrix multiplication even if matrices
#          are not of the same dimension, but are row or
#          column compatible
#---------------------------------------------------
# USAGE: result = matmul(x,y)
# where:    x,y = two matrices (not of the same dimension
#                 but are row or column compatible)
#---------------------------------------------------
# RETURNS: result = x * y where x and y are row or column compatible
# --------------------------------------------------
# written by:
# James P. LeSage, Dept of Economics
# Unioutersity of Toledo
# 2801 W. Bancroft St,
# Toledo, OH 43606
# jpl@jpl.econ.utoledo.edu
matmul <- function(x,y){
  rx <- nrow(x); cx <- ncol(x)
  ry <- nrow(y); cy <- ncol(y)
  if( (cx == cy) && (rx == ry) ){
    results <- x * y
  }else if( (cx == cy) && (rx == 1) ){
    results <- y * repmat(x,ry,1) 
  }else if( (cx == cy) && (ry == 1) ){
    results <- x * repmat(y,rx,1)
  }else if( (rx == ry) && (cx == 1) ){
    results <- y * repmat(x,1,cy)
  }else if( (rx == ry) && (cy == 1) ){
    results <- x * repmat(y,1,cx) 
  }else{
    print('matmul: non-conformable in row or column dimension')
    results <- NULL
  }
  return(results)
}