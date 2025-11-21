
source("Code R/Génération des données_Calcul des vraies valeurs.R")
source("Code R/Fonc01_RegLog.R")
source("Code R/Fonc02_IPW.R")

# Chargements des librairies nécessaires

library(simhelpers)
library(dplyr)
library(tibble)
library(knitr)
library(dplyr)
library(kableExtra)

# Analyser les resultats

Tab01$n <- c("1000", "-", "-", "-")

################################################################################
########################## Régression logistique ###############################

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 8, co_inf_para2 = -8)
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

resultats$vrai_param <- rep(0.6507249, nsim)

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
#|1000 |MCSE_bias |0.02418405            | # |Bias |0.002370728           |
#|-    |MCSE_var  |0.001131156           | # |Var  |0.005848682           |
#|-    |MCSE_mse  |0.001025901           | # |MSE  |0.005269434           |
#|-    |%Cov      |1                     |

################################################################################
################################ IPW ###########################################

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 8, co_inf_para2 = -8)
  resultats2[i,] <- IPW(dat)
  
  # DT : Ajout d'une ligne pour suivre l'avancement
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}

##    - statistiques descriptives

summary(resultats2$RRm)
mean(resultats2$RRm)
sd(resultats2$RRm)

##    - biais, variance, moyenne de l'erreur-type, couverture des IC

resultats2$vrai_param <- rep(0.6326075, nsim) # vraie valeur du paramètre

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
Tab01$IPW[1] <- MCSE_biais[3]

### MCSE_var

MCSE_var <- calc_absolute(resultats2, RRm, vrai_param, criteria = "var")
Tab01$IPW[2] <- MCSE_var[3]


### MCSE_MSE 

MCSE_MSE <- calc_absolute(resultats2, RRm, vrai_param, criteria = "mse")
Tab01$IPW[3] <- MCSE_MSE[3]

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance


coverage <- calc_coverage(resultats2, IC_inf, IC_sup, vrai_param)
Tab01$IPW[4] <- coverage[2] 

kable(Tab01)

#|n    |Methode   |IPW         |
#|:----|:---------|:-----------|
#|1000 |MCSE_bias |0.04395892  | # |Bias |0.1432697   |
#|-    |MCSE_var  |0.006014272 | # |Var  |0.01932387  |
#|-    |MCSE_mse  |0.01506683  | # |MSE  |0.0379177   |
#|-    |%Cov      |0.7         |