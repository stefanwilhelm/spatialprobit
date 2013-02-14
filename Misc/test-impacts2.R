library(spatialprobit)

# marginal effects for standard probit
#
# see http://ideas.repec.org/p/ucn/wpaper/201122.html
# and http://www.r-bloggers.com/probitlogit-marginal-effects-in-r/
mfxboot <- function(modform,dist,data,boot=1000,digits=3){  
  x <- glm(modform, family=binomial(link=dist), data)  
  # get marginal effects  
  pdf <- ifelse(dist=="probit",                
    mean(dnorm(predict(x, type = "link"))),                
    mean(dlogis(predict(x, type = "link"))))  
  marginal.effects <- pdf*coef(x)  
  # start bootstrap  
  bootvals <- matrix(rep(NA,boot*length(coef(x))), nrow=boot)  
  set.seed(1111)  
  for(i in 1:boot){    
    samp1 <- data[sample(1:dim(data)[1],replace=T,dim(data)[1]),]    
    x1 <- glm(modform, family=binomial(link=dist),samp1)    
    pdf1 <- ifelse(dist=="probit",                   
      mean(dnorm(predict(x, type = "link"))),                   
      mean(dlogis(predict(x, type = "link"))))    
    bootvals[i,] <- pdf1*coef(x1)  
  }  
  res <- cbind(marginal.effects,
    apply(bootvals,2,sd),
    marginal.effects/apply(bootvals,2,sd))  
  if(names(x$coefficients[1])=="(Intercept)"){    
    res1 <- res[2:nrow(res),]    
    res2 <- matrix(as.numeric(sprintf(paste("%.",paste(digits,"f",sep=""),sep=""),res1)),nrow=dim(res1)[1])         
    rownames(res2) <- rownames(res1)    
  } else {    
    res2 <- matrix(as.numeric(sprintf(paste("%.",paste(digits,"f",sep=""),sep="")),nrow=dim(res)[1]))    
    rownames(res2) <- rownames(res)    
  }  
  colnames(res2) <- c("marginal.effect","standard.error","z.ratio")    
  return(res2)
}

# number of observations
n <- 100

# true parameters
beta <- c(0, 1, -1)
rho <- 0.75

# design matrix with two standard normal variates as "covariates"
X <- cbind(intercept=1, x1=rnorm(n), x2=rnorm(n))

# sparse identity matrix
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)

# number of nearest neighbors in spatial weight matrix W
m <- 6

# spatial weight matrix with m=6 nearest neighbors
W <- sparseMatrix(i=rep(1:n, each=m),
  j=replicate(n, sample(x=1:n, size=m, replace=FALSE)), x=1/m, dims=c(n, n))

# innovations
eps <- rnorm(n=n, mean=0, sd=1)

# generate data from model
S <- I_n - rho * W
z <- solve(qr(S), X %*% beta + eps)
y <- as.vector(z >= 0)  # 0 or 1, FALSE or TRUE


# estimate (non-spatial) probit model and compute marginal effects
# just to check if our values are reasonable
df <- data.frame(cbind(response=y, X))
colnames(df) <- c("response", "intercerpt", "x1", "x2")
probit.fit1 <- glm(y ~ x1 + x2, data=df, family=binomial(link = "probit"))
summary(probit.fit1)

# predict(probit.fit1, type = "link")
mfx1 <- mfxboot(y ~ x1 + x2, data=df, dist="probit", boot=1000)
mfx1

# estimate SAR probit model
set.seed(12345)
sarprobit.fit1 <- sar_probit_mcmc(y, X, W, ndraw=1000, burn.in=50,
  thinning=1, prior=NULL, computeMarginalEffects=TRUE, showProgress=TRUE)
  
summary(sarprobit.fit1)
impacts(sarprobit.fit1)

set.seed(12345)
mfx <- marginal.effects(sarprobit.fit1)

# means
colMeans(mfx$direct)
colMeans(mfx$indirect)
colMeans(mfx$total)

# quantiles
apply(mfx$direct, 2, quantile, prob=c(0.025, 0.975))

# summaries
mfx$summary_direct
mfx$summary_indirect
mfx$summary_total

# 95% confidence interval of direct impacts based on asymptotic normality
rbind(
 lower=colMeans(mfx$direct) - 1.96 * apply(mfx$direct, 2, sd),
 upper=colMeans(mfx$direct) + 1.96 * apply(mfx$direct, 2, sd)
)

# 95% confidence interval based on quantiles
apply(me$direct, 2, quantile, prob=c(0.025, 0.975))

sarprobit.fit2 <- sarprobit(y ~ x1 + x2, W, data=df, ndraw=1000, burn.in=50,
  thinning=1, prior=NULL, computeMarginalEffects=TRUE, showProgress=TRUE)
  




plot(density(me$direct[,1]))
lines(density(me2$direct[,1]), col="blue")
lines(density(sarprobit.fit1$direct[,1]), col="red")
plot(me$direct[,1],sarprobit.fit1$direct[,1])
abline(a=0, b=1, lty=1, col="red")
# problem: both should be the same, wrong:
# tr(W^i) is determined with simulation, so there will be different realisations in sarprobit.fit1
# and marginal.effects.sarprobit()
quantile(me$direct[,1], prob=c(0.025, 0.975)) # 95% confidence interval
quantile(sarprobit.fit1$direct[,1], prob=c(0.025, 0.975))
}