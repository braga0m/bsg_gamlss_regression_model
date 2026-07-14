#########################################################################pacotes
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse","latex2exp", "Cairo", "gamlss","devtools", "GLMsData",
               "patchwork", "gamlss.ggplots")
####
search()

####diretório
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()


####funções dependentes
source("BSG_GAMLSS.r")

####alterações no plot
newpar <-par(mfrow = c(1,1), #disposição dos gráficos
             mar = par("mar") + c(0,1,0,0),
             col.axis = "black",
             col = "darkgreen", 
             col.main = "black",
             col.lab = "black",
             pch = 1,
             #cex = 0.45, 
             cex.lab = 1.2, 
             cex.axis = 1, 
             cex.main = 1.2, 
             bg = "white")

##########################################################################DENTAL
###
data(dental)
attach(dental)

########## Ajuste da BS Clássica
fit_BS <- gamlss(formula = DMFT ~ Sugar + factor(Indus),
                 sigma.formula = ~ factor(Indus),
                 nu.formula = ~ 1,
                 data = dental,
                 nu.fix = T, 
                 i.nu = 0.5,
                 family = BSG(mu.link = "log", sigma.link = "log"),
                 method = RS())

###########selecão de variáveis
#scope <- list(lower = ~ 1, upper = ~ Sugar + factor(Indus))

#fit_BS_mu <- stepGAIC(fit_BS, what = "mu", scope = scope)
#fit_BS_sigma <- stepGAIC(fit_BS, what = "sigma", scope = scope)
#fit_BS_geral <- stepGAIC(fit_BS, scope = scope)




###########
summary(fit_BS)
plot(fit_BS, par = newpar)
wp(fit_BS)
########## Ajuste da BS de Owen
fit_BSG <- gamlss(formula = DMFT ~ Sugar + factor(Indus),
                  sigma.formula = ~  factor(Indus),
                  nu.formula = ~ 1,
                  data = dental,
                  family = BSG(mu.link = "log"),
                  method = RS())

#########
###########selecão de variáveis
#scope <- list(lower = ~ 1, upper = ~ Sugar + factor(Indus))
#fit_BSG_mu <- stepGAIC(fit_BSG, what = "mu", scope = scope)
#fit_BSG_sigma <- stepGAIC(fit_BSG, what = "sigma", scope = scope)
#fit_BSG_nu <- stepGAIC(fit_BSG, what = "nu", scope = scope)


#########
summary(fit_BSG)
plot(fit_BSG, par = newpar)
wp(fit_BSG)



##########################################################################CHEESE
###
data("cheese")
attach(cheese)

#
#with(cheese,cor(cbind(Taste, Acetic, logH2S=log(H2S), Lactic)))


########## Ajuste da BS Clássica
fit_BS <- gamlss(formula = Taste ~ log(H2S) + Lactic, 
                 sigma.formula = ~ 1,
                 nu.formula = ~ 1,
                 data = cheese, 
                 nu.fix = T, i.nu = 0.5, 
                 family = BSG(mu.link = "log"),
                 method = RS())


###########selecão de variáveis
scope <- list(lower = ~ 1, upper = ~ log(H2S) + Lactic + Acetic)
fit_BS_mu <- stepGAIC(fit_BS, what = "mu", scope = scope)
fit_BS_sigma <- stepGAIC(fit_BS, what = "sigma", scope = scope)
fit_BS_geral <- stepGAIC(fit_BS, scope = scope)


###########
summary(fit_BS)
plot(fit_BS, par = newpar)
wp(fit_BS)

########## Ajuste da BS de Owen
fit_BSG <- gamlss(formula =  Taste ~  log(H2S) + Lactic, 
                  sigma.formula = ~ 1, 
                  nu.formula = ~ 1,
                  data = cheese, 
                  family = BSG(mu.link = "log", sigma.link ="log", nu.link = "logit"),
                  method = RS())


###
###########selecão de variáveis
# scope <- list(lower = ~ 1, upper = ~ log(H2S) + Lactic + Acetic)
# fit_BSG_mu <- stepGAIC(fit_BSG, what = "mu", scope = scope)
# fit_BSG_sigma <- stepGAIC(fit_BSG, what = "sigma", scope = scope)
# fit_BSG_nu <- stepGAIC(fit_BSG, what = "nu", scope = scope)
# fit_BSG_geral <- stepGAIC(fit_BSG, scope = scope)


####
summary(fit_BSG)
plot(fit_BSG, par = newpar)
wp(fit_BSG)



####
cheese2 <- cheese %>% mutate(H2S = log(H2S))
tab <- cor(cheese2)


xtable(tab)
