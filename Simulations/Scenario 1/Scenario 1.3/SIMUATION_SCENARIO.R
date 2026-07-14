####pacotes
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse","data.table","latex2exp", "Cairo", "gamlss", 
               "devtools", "remotes", "RBS")

search()

####pacote RBS e dependência gamlss.nl não disponível no CRAN
#remotes::install_version("gamlss.nl", version = "4.1-0", repos = "http://cran.us.r-project.org")
#devtools::install_github("santosneto/RBS")
#library(RBS)

####diretório
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()


####leitura dos demais scripts
source("BSG_GAMLSS.r")
source("SAMPLE.r")


####nome dos arquivos txt
#
sink(str_c("Results_Scenario", str_sub(getwd(), -3, -1), ".txt"))

#
names_results1 <- c(str_c("Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n40"), 
                    str_c("Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n80"), 
                    str_c("Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n160"),
                    str_c("Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n320"))

#
names_results2 <- c(str_c("SE_Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n40"), 
                    str_c("SE_Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n80"), 
                    str_c("SE_Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n160"),
                    str_c("SE_Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n320"))


names_results3 <- c(str_c("BIAS_Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n40"), 
                    str_c("BIAS_Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n80"), 
                    str_c("BIAS_Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n160"),
                    str_c("BIAS_Estimates_Scenario_", str_sub(getwd(), -3, -1), "_n320"))


#-------------------------------Global Variables-------------------------------#
#definindo os valores dos coef. de regressão 
VN <- c(40, 80, 160, 320) #Sample sizes
VBETA <- c(1.8, -2) #true beta values - MU
VGAMA <- c(-1) #true gamma values - SIGMA
VLAMBDA <- c(log(4)) #true lambda values - NU
VTHETA <- c(VBETA, VGAMA, VLAMBDA) # true theta values
NREP <- 5000 #number of Monte Carlo replicates - 1000 ou 5000 replicas.

#tamanho do vetor de coef de regressão de cada parâmetro. 
kk1 <- length(VBETA) #mu - mediana
kk2 <- length(VGAMA) #sigma - alpha
kk3 <- length(VLAMBDA) #nu - kappa




#---------------------------generated covariates-------------------------------#
#necessário 6 coeficientes de regressão
se1 = 2 ; se2 = 3 #random seed
set.seed(c(se1,se2), kind="Marsaglia-Multicarry") #To ensure repeatability of the experiment
x1 <- runif(VN[1], min = 0, max = 1)#variável explicativa - uniforme U(0,1)
#x2 <- rbinom(VN[1], size = 1, prob = 0.7)#variável explicativa - Bernoulli(0.7)
#x2 <- factor(x2)
#x3 <- rnorm(VN[1], mean = 5, sd = 2)
#x4 <- rpois(VN[1], lambda = 3)#variável explicativa - Possion(3)
X <- matrix(c(rep(1,VN[1]), x1), ncol=2, byrow=F) #regressor matrix for the median submodel
Z <- matrix(c(rep(1,VN[1])), ncol = 1, byrow=F) #regressor matrix for the sigma submodel
W <- matrix(c(rep(1,VN[1])), ncol=1, byrow=F) #regressor matrix for the nu submodel

#dados <- tibble()


#------------------------------Início do loop----------------------------------#
#INÍCIO DO LOOP
for(l in 1:length(VN)){
  N = VN[l]
  
  
  #----------Initializing vectors to store the estimates and statistic z---------#
  thetahat <- matrix(NA, NREP, kk1 + kk2 + kk3)
  se_thetahat <- matrix(NA, NREP, kk1 + kk2 + kk3)
  z_thetahat <- matrix(NA, NREP, kk1 + kk2 + kk3)
  relative_bias <- matrix(NA, NREP, kk1 + kk2 + kk3)
  
  
  #contadores
  cont <- 0 # Counter for Monte Carlo replicates 
  n_descartadas <- 0 # número de replicas descartadas
  
  #
  while(cont < NREP)
  {
    cont <- cont + 1
    perc <- cont/NREP
    
    
    #--------------------To print the progress of the simulation-------------------#
    if(perc == 0.25 || perc == 0.5 || perc == 0.75 || perc ==1) cat("Perc. Replic. MC =",perc*100,"%","\n")
    
    #geração da amostra
    SampleG <- FunSample(n = N, mXini= X, mZini = Z, mWini = W, size = VN[1], theta = c(VBETA, VGAMA, VLAMBDA), linkmu="log",
                         linksigma="log", linknu="logit")
    
    #
    y <- SampleG$y
    mu <- SampleG$mu
    sigma <- SampleG$sigma
    nu <- SampleG$nu
    
    #ajuste para definir um chute inicial mais apropriado
    #alterar de acordo com o cenário
    # fitSTART <- gamlss(formula = y ~ SampleG$X[,2],
    #                    sigma.formula = ~1,
    #                    family = RBS(mu.link = "log", sigma.link = "log"),
    #                    method = CG(),
    #                    trace = FALSE
    #                    )
    
    # fitSTART <- gamlss(formula = y ~ SampleG$X[,2],
    #                    sigma.formula = ~1,
    #                    family = BSG(mu.link = "log", sigma.link = "log", nu.link = "logit"),
    #                    method = RS(),
    #                    trace = FALSE)
    
    #chute inicial
    #mu.inicial <- fitted(fitSTART, what = "mu")
    #sigma.inicial <- fitted(fitSTART, what = "sigma")
    # sigma.inicial <- if_else(sigma.inicial > 10, 1.5, sigma.inicial)
    
    #impedir que o loop pare quando ocorre um erro
    fit_result <- tryCatch({
      n_descartadas <- n_descartadas + 1
      
      #estrutura do modelo
      #alterar a depender do cenário
      fitMLE <- gamlss(formula = y ~ SampleG$X[,2], 
                       sigma.formula = ~ 1, 
                       nu.formula = ~ 1, 
                       family = BSG(), 
                       method = RS(),
                       #mu.start = mu.inicial,
                       #sigma.start = sigma.inicial,
                       #nu.start = rep(0.5, length(y)),
                       trace = FALSE)
      
      #aviso de não convergência
      if (fitMLE$converged && all(!is.na(sqrt(diag(vcov(fitMLE)))))) {
        list(success = TRUE, fit = fitMLE, se = sqrt(diag(vcov(fitMLE))), convergence = fitMLE$converged)
        
        
      }else{
        message("Erro na iteração ", cont, ": o algoritmo não convergiu.")
        list(success = FALSE, convergence = FALSE)
      }
      
      
    }, error = function(e) {
      
      # se erro, imprime mensagem e retorna NULL
      message("Erro na iteração ", cont, ": ", conditionMessage(e))
      list(success = FALSE, fit = FALSE, se = "FALSE", convergence = FALSE)
    })
    
    
    
    #-----------------------parameter estimates------------------------------#
    #
    if (fit_result$success & fit_result$convergence) {
      fitMLE <- fit_result$fit
      std_errors <- fit_result$se
      
      #COEFICIENTES ESTIMADOS DE REGRESSÃO
      thetahat[cont,] <- c(as.numeric(fitMLE$mu.coefficients), 
                           as.numeric(fitMLE$sigma.coefficients), 
                           as.numeric(fitMLE$nu.coefficients))
      
      se_thetahat[cont,] <- c(as.numeric(std_errors[1:ncol(X)]),
                              as.numeric(std_errors[(ncol(X) + 1):(ncol(X) + ncol(Z))]),
                              as.numeric(std_errors[(ncol(X) + ncol(Z) + 1):(ncol(X) + ncol(Z) + ncol(W))]))
      
      z_thetahat[cont, ] <- (thetahat[cont,] - VTHETA) / se_thetahat[cont,]
      
      relative_bias[cont,] <- (thetahat[cont,] - VTHETA)/VTHETA
      
    } else {
      cont <- cont - 1  # descarta e repete essa replicação
      next
    }  
    
    
  }#closes MC replicates
  
  
  
  
  
  
  #-------------------------------RESULTADOS------------------------------------#
  
  
  
  #-----------------------mean of estimates-----------------------#
  M_thetahat <- apply(thetahat,2, function(x) mean(x, na.rm = TRUE))
  
  
  #---------------------median of estimates-----------------------#
  Med_thetahat <- apply(thetahat,2,function(x) median(x, na.rm = TRUE))
  
  #------------------relative bias of estimates-------------------#
  RB_thetahat <- (M_thetahat - VTHETA)/VTHETA
  
  #--------------relative median bias of estimates----------------#
  RMB_thetahat <- (Med_thetahat - VTHETA)/VTHETA
  
  #--------------------mean of s.e. estimates---------------------#
  M_se_thetahat <- apply(se_thetahat,2, function(x) mean(x, na.rm = TRUE))
  
  #--------------------median of s.e. estimate--------------------#
  Med_se_thetahat <- apply(se_thetahat,2,function(x) median(x, na.rm = TRUE))
  
  #--------empirical null levels of Wald-type tests---------------#
  LEVEL_MLE <- apply(abs(z_thetahat) > 1.959964, 2, function(x) mean(x, na.rm = TRUE))
  
  
  
  
  #---------------------------------Outputs--------------------------------------#
  #
  out_estimates_MLE <- cbind(M_thetahat, 
                             RB_thetahat,
                             M_se_thetahat, 
                             LEVEL_MLE, 
                             Med_thetahat, 
                             RMB_thetahat, 
                             Med_se_thetahat)
  
  cat(" ========================================================================================== \n")
  cat(" TRUE VALUES ","\n")
  cat(" Sample size ",N ,"\n")
  
  #
  out_VTHETA <- rbind(VTHETA)
  colnames(out_VTHETA) <- c(str_c("beta", seq(1, (ncol(X)))),
                            str_c("gama", seq(1, (ncol(Z)))),
                            str_c("lambda", seq(1, (ncol(W)))))
  
  
  
  rownames(out_VTHETA) <- c("")
  print(out_VTHETA)
  
  #
  Summ_mu <- round(c(min(mu, na.rm = T), 
                     max(mu, na.rm = T), 
                     median(mu, na.rm = T), 
                     min(sigma, na.rm = T), 
                     max(sigma, na.rm = T), 
                     median(sigma, na.rm = T), 
                     max(sigma)/min(sigma),
                     min(nu, na.rm = T), 
                     max(nu, na.rm = T), 
                     median(nu, na.rm = T),
                     n_desc = n_descartadas - NREP,
                     perc_desc = (n_descartadas - NREP)/NREP),3)
  
  #
  out_Smu <- rbind(Summ_mu)
  colnames(out_Smu)=c("Mu_min", 
                      "Mu_max", 
                      "Mu_median", 
                      "Sigma_min", 
                      "Sigma_max", 
                      "Sigma_median", 
                      "Het_Intensity",
                      "Nu_min", 
                      "Nu_max", 
                      "Nu_median",
                      "n_descartados",
                      "perc_descartados")
  
  rownames(out_Smu)=c("")
  print(out_Smu)
  
  #
  cat(" ========================================================================================== \n")
  cat(" ====================================== MLE estimates ===================================== \n")
  
  
  output_mle = rbind(out_estimates_MLE)
  colnames(output_mle)=c("Mean estimates","Bias","Mean SE", "Null level", "Median estimates","Median Bias","Median SE")
  rownames(output_mle)=c(str_c("MLE_beta", seq(1, (ncol(X)))), 
                         str_c("MLE_gama", seq(1, (ncol(Z)))), 
                         str_c("MLE_lambda", seq(1, (ncol(W)))))
  print(output_mle)
  
  cat(" ========================================================================================== \n")
  
  #------------------------------Print replicates--------------------------------#
  write.table(cbind(thetahat), 
              file = paste(names_results1[l], ".txt", sep = ""), append = FALSE, sep = " ", dec = ".", row.names = FALSE,
              col.names = c(str_c("MLE_beta", seq(1, (ncol(X)))), 
                            str_c("MLE_gama", seq(1, (ncol(Z)))), 
                            str_c("MLE_lambda", seq(1, (ncol(W))))))
  
  write.table(cbind(se_thetahat), 
              file= paste(names_results2[l], ".txt", sep = ""), append = FALSE, sep = " ", dec = ".", row.names = FALSE,
              col.names = colnames(out_VTHETA) <- c(str_c("se_MLE_beta", seq(1, (ncol(X)))), 
                                                    str_c("se_MLE_gama", seq(1, (ncol(Z)))), 
                                                    str_c("se_MLE_lambda", seq(1, (ncol(W))))))
  
  write.table(cbind(relative_bias), 
              file = paste(names_results3[l], ".txt", sep = ""), append = FALSE, sep = " ", dec = ".", row.names = FALSE,
              col.names = c(str_c("bias_MLE_beta", seq(1, (ncol(X)))), 
                            str_c("bias_MLE_gama", seq(1, (ncol(Z)))), 
                            str_c("bias_MLE_lambda", seq(1, (ncol(W))))))
  
  
  
}#closes samples loop

sink()
