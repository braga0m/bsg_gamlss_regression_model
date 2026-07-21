#*************************************************************************Description*************************************************************************************************#
#Replicate X, Z and W of the smallest sample and generate the response variable y.
#***ARGUMENTS***#
# n - sample size. Must be 1, 2, 3, 4, 5, 6, 7, 8, 12, 16 or 32 times the smallest sample size.
# X -  regressor matrix for the mu (median) submodel with number of rows equal to the smallest sample size.
# Z -  regressor matrix for the sigma submodel with number of rows equal to the smallest sample size.
# W -  regressor matrix for the nu submodel with number of rows equal to the smallest sample size.
# size - smallest sample size.
# theta - vector of true parameters.
# linkmu - character specification of the link function in the mu (median) submodel. Currently, "log" and "identity" are supported. Default is "log".
# linksigma - character specification of the link function in the sigma submodel. Currently, "identity", "log", "sqrt" are supported. Default is "log".
# linknu - character specification of the link function in the nu submodel. Currently, "logit" and "probit" are supported. Default is "logit".


#**************************************************************************************************************************************************************************************#

FunSample <- function(n, mXini, mZini, mWini, size, theta, linkmu, linksigma, linknu){
  if(n == size){     
    mX <- mXini
    mZ <- mZini
    mW <- mWini
  }
  else if (n == size*2)
  {
    mX <- rbind(mXini, mXini)
    mZ <- rbind(mZini, mZini)
    mW <- rbind(mWini, mWini)
  }
  else if (n == size*3)
  {
    mX <- rbind(mXini, mXini, mXini)
    mZ <- rbind(mZini, mZini, mZini)
    mW <- rbind(mWini, mWini, mWini)
  }
  else if (n == size*4)
  {
    mX <- rbind(mXini, mXini, mXini, mXini)
    mZ <- rbind(mZini, mZini, mZini, mZini)
    mW <- rbind(mWini, mWini, mWini, mWini)
  }	
  else if (n == size*5)
  {
    mX <- rbind(mXini, mXini, mXini, mXini, mXini)
    mZ <- rbind(mZini, mZini, mZini, mZini, mZini)
    mW <- rbind(mWini, mWini, mWini, mWini, mWini)
  }	
  else if (n == size*6)
  {
    mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini)
    mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini)
    mW <- rbind(mWini, mWini, mWini, mWini, mWini, mWini)
  }
  else if (n == size*7)
  {
    mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini, mXini)
    mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini, mZini)
    mW <- rbind(mWini, mWini, mWini, mWini, mWini, mWini, mWini)
  }
  else if (n == size*8)
  {
    mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini)
    mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini)
    mW <- rbind(mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini)
  }
  else if (n == size*12)
  {
    mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini)
    mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini)
    mW <- rbind(mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini)
  }
  else if (n == size*16)
  {
    mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini)
    mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini)
    mW <- rbind(mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini)
  }
  else if (n == size*32)
  {
    mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini)
    mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini)
    mW <- rbind(mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini, mWini)
  }	
  
  X <- mX; Z <- mZ; W <- mW;  kk1 <- ncol(X); kk2 <- ncol(Z); kk3 <- ncol(W)
  
  # vector of coefficient regressors
  beta <- theta[1:kk1]
  gama <- theta[(kk1+1.0):(kk1+kk2)]
  lambda <- theta[(kk1+kk2+1.0):(kk1+kk2+kk3)]
  
  #linear predictors                                                
  eta <- as.vector(X%*%beta)	
  delta <- as.vector(Z%*%gama)
  zeta <-  as.vector(W%*%lambda) 
  
  if(linkmu == "log") mu <- exp(eta) #log
  if(linkmu == "identity") mu <- eta #identity
  if(linksigma == "sqrt") mu <- eta^2 #sqrt
  if(linksigma == "inverse") mu <- 1/eta #inverse
  
  if(linksigma == "log") sigma <- exp(delta) #log
  if(linksigma == "identity") sigma <- delta #identity
  if(linksigma == "sqrt") sigma <- delta^2 #sqrt
  if(linksigma == "inverse") sigma <-  1/delta #inverse
  
  
  if(linknu == "logit") nu <- exp(zeta)/(1.0+exp(zeta)) #logit
  if(linknu == "probit") nu <- pnorm(zeta) #probit
  if(linknu == "identity") nu <- zeta #identity
  
  #Generating y
  y <- numeric(0)
  
  for(i in 1:n){
    y[i] <- rBSG(1, mu[i], sigma[i], nu[i])                                            
  }
  
  results <- list(y=y, X=X, Z=Z, W=W, mu=mu, sigma= sigma, nu=nu)
  return(results)
}