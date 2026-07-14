#########################################################################pacotes
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse","latex2exp", "Cairo", "gamlss","devtools", "GLMsData",
               "patchwork", "knitr", "kableExtra", "xtable")
search()

####diretório
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

#####
source("BSG_GAMLSS.r")

#####
tema <- theme_bw()+
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text=element_text(size=15),
        plot.title = element_text(hjust = 0.5, size=18))

##########################################################################DENTAL
###########dados
data(dental)
attach(dental)

#####tabela resumo
tab_resumo_dental <- tibble("Mínimo" = min(DMFT),
                            "1° Quartil" = quantile(DMFT, probs = 0.25),
                            "Mediana" = median(DMFT),
                            "Média" = mean(DMFT),
                            "3° Quartil" = quantile(DMFT, probs = 0.75),
                            "Máximo" = max(DMFT)) %>% 
  pivot_longer(cols = everything(),
               names_to = "Medidas", 
               values_to = "DMFT")

#####
xtable(tab_resumo_dental)


#####scatter plot
dental %>% 
  ggplot(aes(x = Sugar, y = DMFT))+
  geom_point(size = 2)+
  tema

#####box-plot
dental %>% 
  ggplot(aes(x = Indus, y = DMFT))+
  geom_boxplot(fill = "#def2ff")+
  labs(x = "")+
  tema

######box-plot y
dental %>% 
  ggplot(aes(x = "", y = DMFT))+
  geom_boxplot(fill = "#def2ff")+
  labs(x = "")+
  tema


##########################################################################CHEESE
###########
data("cheese")
attach(cheese)



#####tabela resumo
tab_resumo_cheese <- tibble("Mínimo" = min(Taste),
                            "1° Quartil" = quantile(Taste, probs = 0.25),
                            "Mediana" = median(Taste),
                            "Média" = mean(Taste),
                            "3° Quartil" = quantile(Taste, probs = 0.75),
                            "Máximo" = max(Taste)) %>% 
  pivot_longer(cols = everything(),
               names_to = "Medidas", 
               values_to = "Valores")


#####scatter plot
cheese %>% 
  ggplot(aes(x = Acetic, y = Taste))+
  geom_point(size = 2)+
  tema

cheese %>% 
  ggplot(aes(x = log(H2S), y = Taste))+
  geom_point(size = 2)+
  tema

cheese %>% 
  ggplot(aes(x = Lactic, y = Taste))+
  geom_point(size = 2)+
  tema


######box-plot y
cheese %>% 
  ggplot(aes(x = "", y = Taste))+
  geom_boxplot(fill = "#def2ff")+
  labs(x = "")+
  tema


