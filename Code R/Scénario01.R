
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
  
  vraiRRc <- glm(Y ~ V + C, family = binomial(link = "logit"), data = dat_complet)
  
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
# Initialiser des objets pour contenir les resultats

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))


################################################################################
########################## Analyse des résultats ###############################

nsim <- 100

# Initialiser des objets pour contenir les resultats

resultats <- data.frame(matrix(ncol = 4, 
                               nrow = nsim))

colnames(resultats) <- c("coe_reg", 
                         "err_reg",
                         "RRc", # Risque relatif conditionnel
                         "est_VE")

resultats2 <- data.frame(matrix(ncol = 5, 
                                nrow = nsim))

colnames(resultats2) <- c("RRm",# Risque relatif marginal
                          "VE",
                          "var_log_RRm",# Variance du log du risque relatif marginal
                          "IC_inf", # Borne inférieure de l'intervalle de confiance
                          "IC_sup") # Sa borne supérieure


for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 8, co_inf_para2 = -8, popsize = 1*10**6)
  
  resultats[i,] <- RegLog(dat) # Régression logistique
  resultats2[i,] <- IPW(dat)   # IPW

  # DT : Ajout d'une ligne pour suivre l'avancement
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}
########################## Régression logistique ###############################

Tab01$Methode <- c("RegLog", "-", "-","-")

##    - statistiques descriptives

summary(resultats$RRc)
mean(resultats$RRc)
sd(resultats$RRc)
sd(log(resultats$RRc));
mean(resultats$err_reg);

help("calc_absolute") # calculer les différentes mesures de performance

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats$vrai_param <- rep(mean(l_vraiRRc), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats, RRc, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats, RRc, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats, RRc, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

resultats$lim_inf <- exp(resultats[, 1] - 1.96*resultats[, 2]) 
resultats$lim_sup <- exp(resultats[, 1] + 1.96*resultats[, 2])

mean(resultats$lim_inf < resultats$vrai_param & resultats$lim_sup > resultats$vrai_param)

coverage <- calc_coverage(resultats, lim_inf, lim_sup, vrai_param)

Tab01$`Erreur de Monte Carlo` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais[3]),
  list(name = "MCSE_var", value = MCSE_var[3]),
  list(name = "MCSE_mse", value = MCSE_mse[3]),
  list(name = "%Cov", value = coverage[2])
  
)

Tab01$`Autres` = list(
  
  list(name = "Bias", value = MCSE_biais[2]),
  list(name = "Var", value = MCSE_var[3]),
  list(name = "Mse", value = MCSE_mse[3]),
  list(name = "Précision_var", value = sd(log(resultats$RRc)) - mean(resultats$err_reg)) # Voir si la variance est bien estimée
  
)

kable(Tab01)

#|n    |Methode |Erreur de Monte Carlo                      |Autres                                     |
#|:----|:-------|:------------------------------------------|:------------------------------------------|
#|1000 |RegLog  |MCSE_bias         , 0.0071268926221513     |Bias               , 0.00428675767604619   |
#|-    |-       |MCSE_var            , 0.000664559421799917 |Var                 , 0.000664559421799917 |
#|-    |-       |MCSE_mse            , 0.000685994345421483 |Mse                 , 0.000685994345421483 |
#|-    |-       |%Cov, 0.96                                 |Précision_var     , 0.0186125857979366     |

# K_coverage coverage coverage_mcse width width_mcse
# <int>    <dbl>         <dbl> <dbl>      <dbl>
#   1       1000    0.936       0.00774 0.429    0.00148
# > mean(resultats$err_reg)
# > sd(log(resultats$RRc))

################################ IPW ###########################################

Tab01$Methode <- c("IPW", "-", "-","-")

##    - statistiques descriptives

summary(resultats2$RRm)
mean(resultats2$RRm)
sd(resultats2$RRm)

##    - biais, variance, moyenne de l'erreur-type, couverture des IC

resultats2$vrai_param <- rep(mean(l_vraiRRm), nsim) # vraie valeur du paramètre

Tab01$Methode <- c("IPW")

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats2, RRm, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats2, RRm, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_MSE <- calc_absolute(resultats2, RRm, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance


coverage <- calc_coverage(resultats2, IC_inf, IC_sup, vrai_param)

Tab01$`Erreur de Monte Carlo` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais[3]),
  list(name = "MCSE_var", value = MCSE_var[3]),
  list(name = "MCSE_mse", value = MCSE_mse[3]),
  list(name = "%Cov", value = coverage[2])
  
)

Tab01$`Autres` = list(
  
  list(name = "Bias", value = MCSE_biais[2]),
  list(name = "Var", value = MCSE_var[3]),
  list(name = "Mse", value = MCSE_mse[3]),
  list(name = "Précision_var", value = sd(log(resultats2$RRm)) - mean(sqrt(resultats2$var_log_RRm)) # Voir si la variance est bien estimée
       
  )
  
)
  


kable(Tab01)

#|n    |Methode |Erreur de Monte Carlo                      |Autres                                     |
#|:----|:-------|:------------------------------------------|:------------------------------------------|
#|1000 |IPW     |MCSE_bias          , 0.00714109026683143   |Bias               , 0.00366123637127513   |
#|-    |IPW     |MCSE_var            , 0.000628107910234781 |Var                 , 0.000628107910234781 |
#|-    |IPW     |MCSE_mse            , 0.000685994345421483 |Mse                 , 0.000685994345421483 |
#|-    |IPW     |%Cov, 0.93                                 |Précision_var     , 0.0298310159785812     |

# K_coverage coverage coverage_mcse width width_mcse
# <int>    <dbl>         <dbl> <dbl>      <dbl>
#   1       1000    0.936       0.00774 0.429    0.00148
# > mean(sqrt(resultats2$var_log_RRm))
# [1] 0.1554226
# > sd(log(resultats2$RRm))
# [1] 0.1662523