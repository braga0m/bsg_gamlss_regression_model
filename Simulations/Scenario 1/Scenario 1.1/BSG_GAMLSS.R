####pacotes
if (!require("pacman")) install.packages("pacman")
pacman::p_load("\mutidyverse","data.table","latex2exp", "Cairo", "gamlss")




#-------------------------------------------------------------------------------
# #expressão da função densidade de probabilidade, necessária para obter numericamente as derivadas
# l_BSG <- expression(log(
#   (1/sqrt(2*pi))*(1/(sigma*sqrt(mu)*(y^nu)))*(1-nu+((nu*mu)/y))*exp(-(1/(2*sigma^2))*(((y-mu)^2)/(mu*y^(2*nu))))
# ))
# 
# ######derivadas em relação a cada parâmetro
# #derivadas de primeira ordem
# m1 <- D(l_BSG, "mu")
# s1 <- D(l_BSG, "sigma")
# n1 <- D(l_BSG, "nu")
# 
# 
# 
# #derivadas de segunda ordem
# ms2 <- D(m1, "sigma") #igual a sm2
# mn2 <- D(m1, "nu") #igual a nm2
# sn2 <- D(s1, "nu") #igual ns2


#------------------------------------------------------------------------------
# Utilizando a distribuicao BSG de Owen no GAMLSS
#

#mu = mu = mediana
#sigma = alpha
#nu = kappa = relacionado com a correlação dos ciclos
BSG <- function(mu.link = "log", sigma.link = "log", nu.link = "logit"){
  
  mstats <- checklink("mu.link", "BSG", substitute(mu.link),
                      c("inverse", "log", "identity", "own", "sqrt"))
  
  dstats <- checklink("sigma.link", "BSG", substitute(sigma.link),
                      c("inverse", "log", "identity", "own", "sqrt"))
  
  vstats <- checklink("nu.link", "BSG", substitute(nu.link),
                      c("logit", "probit", "cloglog", "cauchit", "own", "identity"))
  
  
  structure(
    list(family = c("BSG", "Birnbaum-Saunders Generalizada"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE),
         nopar = 3, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         nu.link = as.character(substitute(nu.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         nu.linkfun = vstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         nu.dr = vstats$mu.eta,
         
         
         ######derivadas de primeira ordem
         dldm =  function(y,mu, sigma, nu){
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "mu", delta = 1e-04), "gradient"))
           dldm
         },
         
         dldd = function(y,mu, sigma, nu){
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           dldd
         },
         
         dldv = function(y,mu, sigma, nu){
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "nu", delta = 1e-04), "gradient"))
           dldv
         },
         
         
         #####derivadas de segunda ordem - mesmo parâmetro
         d2ldm2 = function(y,mu, sigma, nu) {
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "mu", delta = 1e-04), "gradient"))
           #sapply(seq_along(y), function(i) {
           #hess <- hessian(
           #func = function(par) dBGEV(y[i], mu = par, sigma = sigma[i], nu = nu[i], tau = tau[i], log = TRUE),
           #x = mu[i]
           #)
           #hess[1, 1]  # segunda derivada em relação a mu
           #})
           d2ldm2 = -dldm^2
           d2ldm2
         }, 
         
         d2ldd2 = function(y,mu, sigma, nu) {
             dldd = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "sigma", delta = 1e-04), "gradient"))
             d2ldd2 = -dldd^2
             d2ldd2
         },
         
         d2ldv2 = function(y,mu, sigma, nu) {
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "nu", delta = 1e-04), "gradient"))
           d2ldv2 = -dldv^2
           d2ldv2
         }, 
         
         #####derivadas de segunda ordem - parâmetro cruzados
         d2ldmdd = function(y,mu, sigma, nu) {
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "mu", delta = 1e-04), "gradient"))
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           
           d2ldmdd = -(dldm*dldd)
           #d2ldmdd = ifelse(d2ldmdd < -1e-15, d2ldmdd,-1e-15)#add
           d2ldmdd
         },
         
         d2ldmdv = function(y,mu, sigma, nu) {
           dldm = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "mu", delta = 1e-04), "gradient"))
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "nu", delta = 1e-04), "gradient"))
           
           d2ldmdv = -(dldm*dldv)
           #d2ldmdv = ifelse(d2ldmdv < -1e-15, d2ldmdv,-1e-15)#add
           d2ldmdv
         },
         
         d2ldddv = function(y,mu, sigma, nu) {
           dldd = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "sigma", delta = 1e-04), "gradient"))
           dldv = as.vector(attr(gamlss:::numeric.deriv(dBSG(y, mu, sigma, nu, log = TRUE), "nu", delta = 1e-04), "gradient"))
           d2ldddv = -(dldd*dldv)
           #d2ldddv = ifelse(d2ldddv <-1e-15, d2ldddv,-1e-15)#add
           d2ldddv
         },
         
         
         #####
         G.dev.incr = function(y,mu,sigma,nu, ...) -2*dBSG(y = y, mu = mu, sigma = sigma , nu = nu, log=TRUE),
         rqres = expression(rqres(pfun="pBSG", type = "Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
         
         #####valores inicias para mu e alpha
         mu.initial = expression(mu <-  rep(median(y),length(y))),
         sigma.initial = expression(sigma <- rep(1,length(y))),
         nu.initial = expression(nu <- rep(0.5, length(y))),
         
         #####restrição ao domínio dos parâmetros e da v.a
         mu.valid = function(mu) all(is.finite(mu) &  mu > 0) ,
         #sigma.valid = function(sigma) all(is.finite(sigma) & (sigma > 0 & sigma < 2)),
         sigma.valid = function(sigma) all(is.finite(sigma) & (sigma > 0)),
         nu.valid = function(nu) all(is.finite(nu) & (nu > 0.01 & nu < 0.99)),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}








#-------------------------------------------------------------------------------
# fdp da distribuicao BSG - OK!
dBSG<-function(y, mu=1, sigma=1, nu = 0.5, log=FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  #if (any(nu <= 0 | nu >= 1)) stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(nu <= 0.01 | nu >= 0.99)) stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(y < 0)) stop(paste("y must be positive", "\n", ""))
  
  fy1 <- (1/sqrt(2*pi))*(1/(sigma*sqrt(mu)*(y^nu)))*(1-nu+((nu*mu)/y))*exp(-(1/(2*sigma^2))*(((y-mu)^2)/(mu*y^(2*nu))))
  
  if(log == TRUE) fy <- log(fy1) else fy <- fy1
  #fy <- ifelse(y <= 0, 0, fy)
  
  return(fy)
}







#-------------------------------------------------------------------------------
# fda da distribuicao BS - OK!
pBSG <- function(q, mu=1, sigma=1, nu=0.5, lower.tail = TRUE, log = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  #if (any(nu <= 0) | any(nu >= 1)) stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(nu <= 0.01 | nu >= 0.99)) stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(q <= 0)) stop(paste("y must be positive", "\n", ""))
  
  #f_arg <- (1/sigma)*((y^(1 - nu))/sqrt(mu) - sqrt(mu)/(y^nu))
  #cdf <- pnorm(f_arg, mean = 0, sd = 1)
  cdf <- pnorm((1/sigma)*((q^(1-nu))/sqrt(mu)-sqrt(mu)/(q^nu)))  
  
  if(lower.tail==TRUE) cdf <- cdf else cdf <- 1 - cdf
  if(log==TRUE) cdf <- log(cdf)
  
  return(cdf)
}

#-------------------------------------------------------------------------------
#função quantílica 
qBSG <- function(p, mu=1, sigma=1, nu = 0.5, lower.tail = TRUE, log = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  #if (any(nu <= 0) | any(nu >= 1)) stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(nu <= 0.01 | nu >= 0.99)) stop(paste("nu must be between 0 and 1", "\n", ""))
  if (log==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  root=NULL
  for(i in 1:length(p)){
    
    prob=p[i]
    
    #f=function(y_aux,sigma,mu,nu,p) pBSG(y_aux,sigma,mu,nu)- p
    f=function(y, mu, sigma, nu, p) pBSG(y, mu, sigma, nu) - p
    
    #root[i]=uniroot(f,c(0.0000001, 10000000), tol = 0.0001, sigma, mu, nu, p)$root
    root[i]=uniroot(f,c(0.0000001, 10000000), tol = 0.0001, mu, sigma, nu, p)$root
  }
  return(root)
}


#-------------------------------------------------------------------------------
#gerador de números aleatórios
rBSG <- function(n, mu = 1, sigma = 1, nu = 0.5){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  #if (any(nu <= 0) | any(nu >= 1)) stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(nu <= 0.01 | nu >= 0.99)) stop(paste("nu must be between 0 and 1", "\n", ""))
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  
  y <- numeric(0)
  
  # Gerando amostras GBS de tamanho n
  for(i in 1:n) 
  {
    z <-  rnorm(1)
    f <- function(y, sigma, mu, nu) sigma*z*(y^nu)*sqrt(mu)+mu-y
    
    y[i] <- uniroot(f,c(0.000001, 1000000000000000), tol = 0.0001, sigma, mu, nu)$root
  }
  
  return(y)
}



