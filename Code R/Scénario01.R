
source("Code R/Génération des données.R")
source("Code R/Fonc01_RegLog.R")
source("Code R/Fonc02_IPW.R")

# Chargements des librairies nécessaires

library(simhelpers)
library(dplyr)
library(tibble)
library(knitr)
library(dplyr)
library(kableExtra)

# Calcul des vraies valeurs 

set.seed(1) # Pour avoir toujours les mêmes germes

nsim <- 10

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000,
                      co_inf_para1 = 8, co_inf_para2 = -8)
  
  summary(dat)
  
  dat0 <- data.frame(C = dat$C, V = 0, Y = dat$I2_0*dat$W2_0*dat$H_0)
  
  dat1 <- data.frame(C = dat$C, V = 1, Y = dat$I2_1*dat$W2_1*dat$H_1) 
  
  dat_complet <- rbind(dat0, dat1)
  
  vraiRRc <- glm(Y ~ V + C, family = binomial(link = "log"), data = dat_complet)
  
  l_vraiRRc[i] <- exp(coef(vraiRRc)[2])
  
  l_vraiRRm[i] <- mean(dat1$Y)/mean(dat0$Y)
  
  print(data.frame(Sys.time(), i))
  
}

l_vraiRRc

mean(l_vraiRRc)
sd(l_vraiRRc)


l_vraiRRm

mean(l_vraiRRm)
sd(l_vraiRRm)

################################################################################
# Analyser les resultats

Tab01 <- data.frame(matrix(ncol = 3, 
                           nrow = 4))

Tab01$n <- c("1000", "-", "-", "-")

colnames(Tab01) <- c("n",
                     "Methode", 
                     "Regression_logistique")

Tab01$Methode <- c("MCSE_bias", "MCSE_var", "MCSE_mse", "%Cov")

################################################################################
########################## Régression logistique ###############################

nsim <- 1000

# Initialiser des objets pour contenir les resultats

resultats <- data.frame(matrix(ncol = 6, 
                               nrow = nsim))

colnames(resultats) <- c("coe_reg", 
                         "err_reg",
                         "RRc", # Risque relatif conditionnel
                         "est_VE")

resultats2 <- data.frame(matrix(ncol = 5, 
                                nrow = nsim))

colnames(resultats2) <- c("RRm",# Risque relatif marginal
                          "VE",
                          "var_RRm",# Variance du risque relatif marginal
                          "IC_inf", # Borne inférieure de l'intervalle de confiance
                          "IC_sup") # Sa borne supérieure


for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 8, co_inf_para2 = -8, popsize = 1*10**6)
  resultats[i,] <- RegLog(dat)

  # DT : Ajout d'une ligne pour suivre l'avancement
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}
##    - statistiques descriptives

summary(resultats$RRc)
mean(resultats$RRc)
sd(resultats$RRc)

help("calc_absolute") # calculer les différentes mesures de performance

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats$vrai_param <- rep(mean(l_vraiRRc), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats, RRc, vrai_param, criteria = "bias")
Tab01$Regression_logistique[1] <- MCSE_biais[3]

### MCSE_var

MCSE_var <- calc_absolute(resultats, RRc, vrai_param, criteria = "var")
Tab01$Regression_logistique[2] <- MCSE_var[3]

### MCSE_MSE 

MCSE_MSE <- calc_absolute(resultats, RRc, vrai_param, criteria = "mse")
Tab01$Regression_logistique[3] <- MCSE_MSE[3]

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

resultats$lim_inf <- exp(resultats[, 1] - 1.96*resultats[, 2]) 
resultats$lim_sup <- exp(resultats[, 1] + 1.96*resultats[, 2])

mean(resultats$lim_inf < resultats$vrai_param & resultats$lim_sup > resultats$vrai_param)

coverage <- calc_coverage(resultats, lim_inf, lim_sup, vrai_param)
Tab01$Regression_logistique[4] <- coverage[2] 

kable(Tab01)

#|n    |Methode   |Regression_logistique | 
#|:----|:---------|:---------------------|
#|1000 |MCSE_bias |0.03777883            | # |Bias |-0.01649214           |
#|-    |MCSE_var  |0.003786253           | # |Var  |0.0142724             |
#|-    |MCSE_mse  |0.004233154           | # |MSE  |0.01311715            |
#|-    |%Cov      |0.9                   |

################################################################################
################################ IPW ###########################################

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 8, co_inf_para2 = -8, popsize = 10**6)
  resultats2[i,] <- IPW(dat)
  
  # DT : Ajout d'une ligne pour suivre l'avancement
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}

##    - statistiques descriptives

summary(resultats2$RRm)
mean(resultats2$RRm)
sd(resultats2$RRm)

##    - biais, variance, moyenne de l'erreur-type, couverture des IC

resultats2$vrai_param <- rep(mean(l_vraiRRm), nsim) # vraie valeur du paramètre

Tab01 <- data.frame(matrix(ncol = 3, 
                           nrow = 4))

colnames(Tab01) <- c("n",
                     "Methode", 
                     "IPW")

Tab01$n <- c("1000", "-", "-", "-")
Tab01$Methode <- c("MCSE_bias", "MCSE_var", "MCSE_mse", "%Cov")

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats2, RRm, vrai_param, criteria = "bias")
Tab01$IPW[1] <- MCSE_biais[2]

### MCSE_var

MCSE_var <- calc_absolute(resultats2, RRm, vrai_param, criteria = "var")
Tab01$IPW[2] <- MCSE_var[2]


### MCSE_MSE 

MCSE_MSE <- calc_absolute(resultats2, RRm, vrai_param, criteria = "mse")
Tab01$IPW[3] <- MCSE_MSE[2]

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance


coverage <- calc_coverage(resultats2, IC_inf, IC_sup, vrai_param)
Tab01$IPW[4] <- coverage[2] 

kable(Tab01)

#|n    |Methode   |IPW         |
#|:----|:---------|:-----------|
#|1000 |MCSE_bias |0.04684841  | # |Bias |0.05237775   |
#|-    |MCSE_var  |0.008997592 | # |Var  |0.02194774   |
#|-    |MCSE_mse  |0.01295461  | # |MSE  |0.02249639   |
#|-    |%Cov      |0.9         |