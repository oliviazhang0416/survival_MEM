require(MASS)
library(randomizr)

#TE: hazard ratio
gen_surv_dat <- function(n, 
                         dta_description, 
                         rho = 0.1, 
                         outcome_type = "surv", 
                         TE, 
                         censor_t, 
                         N0, N1,
                         source,
                         seeds){
  p <- nrow(dta_description)
  S <- dta_description$sd
  mu <- dta_description$mu
  theta <- dta_description$theta
  binary <- dta_description$binary
  b0 <- unique(dta_description$intercept)
  b <- c(dta_description$coef,log(TE))
  
  
  # covariance matrix
  Vmat <- matrix(0, p, p)
  for(i in 1:ncol(Vmat)){
    for(j in 1:nrow(Vmat)){
      if(i==j){Vmat[i, j] <- S[i]^2}
      if(i!=j){Vmat[i, j] <- rho*S[i]*S[j]}
    }
  }
  set.seed(seeds)
  xmat <- mvrnorm(n = n, mu = mu, Sigma = Vmat)
  
  for(i in 1:p){
    if(binary[i]){
      x <- ifelse(xmat[, i]<qnorm(theta[i]), 1, 0)
      xmat[, i] <- x
    }
  }
  
  ####randomize
  if (source == "internal"){
    set.seed(seeds)
    #Z <- c(rbinom(60,size=1,prob=2/3))
    Z <- complete_ra(N = n,  m = N1)
    xmat <- cbind(xmat,Z)
  }
  else{
    Z <- rep(0,n)
    xmat <- cbind(xmat,Z) 
  }
  
  if(outcome_type == "surv"){
    lpd <- b0 + xmat%*%b 
    set.seed(seeds)
    
    y <- rweibull(n = n, shape = 1, scale = exp(lpd))
    d <- ifelse(y <= censor_t, 1, 0)
    t <- ifelse(d==1, y, censor_t)
  }
  dta <- data.frame(cbind(xmat, exp(lpd), t, d))
  names(dta) <- c(dta_description$xname, "Z","exp(lpd)","t", "d")
  return(dta)
}

