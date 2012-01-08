#PURPOSE: print the results of the sar probit estimation via mcmc
# --------------------------------------------------------------------
#INPUT:
#     object     = structure with the results from sar_probit_mcmc
#     var_names  = vector with names for the parameters under analysis
# --------------------------------------------------------------------
#RETURNS: this functions does not return any values
# --------------------------------------------------------------------
summary.sarprobit <- function(object, var_names=NULL, file=NULL, ...){
  # TODO: check for class "sarprobit"
  if (!inherits(object, "sarprobit")) 
        stop("use only with \"sarprobit\" objects")
        
  nobs      <- object$nobs
  nvar      <- object$nvar
  ndraw     <- object$ndraw
  nomit     <- object$nomit
  #bayesian estimation
  bout_mean <- as.matrix(c(apply(object$bdraw,2,mean),mean(object$pdraw))) #parameter mean column
  bout_sd   <- as.matrix(c(apply(object$bdraw,2,sd)  ,sd(object$pdraw))) #parameter sd colum
  bout_sig  <- matrix(data=NA, nrow=nrow(bout_mean),ncol=1)
  #build bayesian significance levels
  draws     <- cbind( object$bdraw, object$pdraw )
  for( i in 1:ncol(draws) ){
    if( bout_mean[i,1] > 0){
      cnt <- which( draws[,i] > 0 )
    }else{
      cnt <- which( draws[,i] < 0 )
    }
    bout_sig[i,1] <- 1 - (length(cnt)/(ndraw-nomit))
  }
  #standar assymptotic measures
  bout_t    <- bout_mean / bout_sd             #t score b/se
  bout_tPval<- (1 - pt( abs(bout_t), nobs ))*2 #two tailed test = zero probability = z-prob
  #name definition
  if( is.null(var_names)){
    bout_names<- as.matrix(object$names)
  }else{
    bout_names<- as.matrix(var_names)
  }
  
  if(is.null(file)){file <- ""}#output to the console
  #HEADER
  write(sprintf("------------MCMC spatial autoregressive probit------------"), file, append=T)
  #sprintf("Dependent Variable")
  write(sprintf("Execution time  = %6.3f", object$time)  , file, append=T)
  write(sprintf("N steps for TMVN= %6d"  , object$nsteps), file, append=T)
  write(sprintf("N draws         = %6d, N omit (burn-in)= %6d", ndraw, nomit), file, append=T)
  write(sprintf("N observations  = %6d, K covariates    = %6d", nobs, nvar)  , file, append=T)
  write(sprintf("# fo 0 Y values = %6d, # of 1 Y values = %6d", object$zip, nobs - object$zip) , file, append=T)
  write(sprintf("Min rho         = % 6.3f, Max rho         = % 6.3f", object$rmin, object$rmax), file, append=T)
  write(sprintf("----------------------------------------------------------"), file, append=T)
  write(sprintf(""), file, append=T)
  #ESTIMATION RESULTS
  write(sprintf("%30s   %15s %15s %15s %15s %15s", 
                'Parameter Name', 
                'Coefficient', 
                'Standard Dev', 
                'Bayes p-level', 
                't-score',
                'z-prob'), file, append=T)
  for( i in 1:nrow(bout_mean)){
    write(sprintf("%30s   % 15.4f % 15.4f % 15.4f % 15.4f % 15.4f", 
                  bout_names[i,1], 
                  bout_mean[i,1], 
                  bout_sd[i,1], 
                  bout_sig[i,1],
                  bout_t[i,1],
                  bout_tPval[i,1]
          ), file, append=T)          
  }
}