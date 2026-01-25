
source("Code R/Génération des données.R")
source("Code R/Fonc01_RegLog.R")
source("Code R/Fonc02_IPW.R")
source("Code R/Fonc03_TNDDR.R")

# Chargements des librairies nécessaires

library(simhelpers)
library(dplyr)
library(tibble)
library(knitr)
library(dplyr)
library(kableExtra)

# Calcul des vraies valeurs 

set.seed(1) # Pour avoir toujours les mêmes germes

nsim <- 500 # Nombre de réplications

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

resultats3 <- data.frame(matrix(ncol = 9, 
                                nrow = nsim))

colnames(resultats3) <- c("RRm",# Risque relatif marginal
                          "VE",
                          "var_log_RRm",# Variance du log du risque relatif marginal
                          "IC_inf1", # Borne inférieure du premier intervalle de confiance
                          "IC_sup1", # Sa borne supérieure
                          "IC_inf2", # Borne inférieure du deuxième intervalle de confiance
                          "IC_sup2", # Sa borne supérieure
                          "IC_inf3", # Borne inférieure du troisième intervalle de confiance
                          "IC_sup3")# Sa borne supérieure

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 8, co_inf_para2 = -8, popsize = 1*10**6)
  
  resultats[i,] <- RegLog(dat) # Régression logistique
  resultats2[i,] <- IPW(dat)   # IPW
  resultats3[i,] <- TNDDR(dat, RandomForest) # Estimateur doublement robuste
  
  # DT : Ajout d'une ligne pour suivre l'avancement
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}


########################## Régression logistique ###############################

Tab01$Methode <- c("RegLog", "-", "-","-")

##    - statistiques descriptives

summary(resultats$RRc)
mean(resultats$RRc)
sd(resultats$RRc)

sd(resultats$coe_reg)
mean(resultats$err_reg)

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
  list(name = "Précision_var", value = sd(resultats$coe_reg) - mean(resultats$err_reg)) # Voir si la variance est bien estimée
  
)

kable(Tab01)

#|n    |Methode |Erreur de Monte Carlo                      |Autres                                     |
#|:----|:-------|:------------------------------------------|:------------------------------------------|
#|1000 |RegLog  |MCSE_bias          , 0.00295418588667984   |Bias               , 0.00729010837813915   |
#|-    |-       |MCSE_var           , 0.000284378344670318  |Var                , 0.000284378344670318  |
#|-    |-       |MCSE_mse           , 0.000298196850968069  |Mse                , 0.000298196850968069  |
#|-    |-       |%Cov               , 0.956                 |Précision_var      , -0.00133713826292561  |

#K_coverage coverage coverage_mcse width width_mcse
#<int>    <dbl>         <dbl> <dbl>      <dbl>
#  1        500    0.956       0.00917 0.265    0.00181
# > mean(resultats$err_reg)
# > sd(resultats$coe_reg))

################################ IPW ###########################################

Tab01$Methode <- c("IPW", "-", "-","-")

##    - statistiques descriptives

summary(resultats2$RRm)
mean(resultats2$RRm)
sd(resultats2$RRm)

mean(sqrt(resultats2$var_log_RRm))
sd(log(resultats2$RRm))

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
  list(name = "Var", value = MCSE_var[2]),
  list(name = "Mse", value = MCSE_mse[2]),
  list(name = "Précision_var", value = sd(log(resultats2$RRm)) - mean(sqrt(resultats2$var_log_RRm)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)

#|n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
#|:----|:-------|:------------------------------------------|:----------------------------------------|
#|1000 |IPW     |MCSE_bias          , 0.00309524943623309   |Bias               , 0.00617727884141961 |
#|-    |IPW     |MCSE_var           , 0.000337155855902644  |Var                , 0.00479028453625063 |
#|-    |IPW     |MCSE_mse           , 0.000298196850968069  |Mse                , 0.00440802559244113 |
#|-    |IPW     |%Cov               , 0.93                  |Précision_var      , 0.0089665270230341  |

#  K_coverage coverage coverage_mcse width width_mcse
#<int>    <dbl>         <dbl> <dbl>      <dbl>
#  1        500     0.93        0.0114 0.258    0.00170
# > mean(sqrt(resultats2$var_log_RRm))
# [1] 0.1492768
# > sd(log(resultats2$RRm))
# [1] 0.1582433

########################## TNDDR ###############################

Tab01$Methode <- c("TNDDR", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm)
mean(resultats3$RRm)
sd(resultats3$RRm)

mean(sqrt(resultats3$var_log_RRm))
sd(log(resultats3$RRm))

help("calc_absolute") # calculer les différentes mesures de performance

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats3$vrai_param <- rep(mean(l_vraiRRc), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, IC_inf2, IC_sup2, vrai_param)

Tab01$`Erreur de Monte Carlo` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais[3]),
  list(name = "MCSE_var", value = MCSE_var[3]),
  list(name = "MCSE_mse", value = MCSE_mse[3]),
  list(name = "%Cov", value = coverage[2])
  
)

Tab01$`Autres` = list(
  
  list(name = "Bias", value = MCSE_biais[2]),
  list(name = "Var", value = MCSE_var[2]),
  list(name = "Mse", value = MCSE_mse[2]),
  list(name = "Précision_var", value = sd(log(resultats3$RRm)) - mean(sqrt(resultats3$var_log_RRm)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)

# Avec la méthode de la forêt aléatoire

#|n    |Methode |Erreur de Monte Carlo                      |Autres                                     |
#|:----|:-------|:------------------------------------------|:------------------------------------------|
#|1000 |TNDDR   |MCSE_bias          , 0.00315932660029385   |Bias               , 0.0160141383783473    |
#|-    |-       |MCSE_var           , 0.000331549017625076  |Var                , 0.00499067228366214   |
#|-    |-       |MCSE_mse           , 0.000372952913256951  |Mse                , 0.00523714356709567   |
#|-    |-       |%Cov               , 0.968                 |Précision_var      , -0.00595712229257123  | # Avec CI 1
#|-    |-       |%Cov               , 0.968                 |Précision_var      , 0.0789107904397003    | # Avec CI 2
#|-    |-       |%Cov               , 0.968                 |Précision_var      , 0.0789107904397003    | # Avec CI 3

# Régression Lasso

#|n    |Methode |Erreur de Monte Carlo                    |Autres                                 |
#|:----|:-------|:----------------------------------------|:--------------------------------------|
#|1000 |TNDDR   |MCSE_bias         , 0.0436181384545875   |Bias              , 0.0430271139223537 |
#|-    |-       |MCSE_var          , 0.00662407462016709  |Var               , 0.0190254200224357 |
#|-    |-       |MCSE_mse          , 0.00951707242118954  |Mse               , 0.0189742105526793 |
#|-    |-       |%Cov              , 0.9                  |Précision_var     , 0.0651750510281717 | # Avec CI 1
#|-    |-       |%Cov              , 0.9                  |Précision_var     , 0.0651750510281717 | # Avec CI 2
#|-    |-       |%Cov              , 0.9                  |Précision_var     , 0.0651750510281717 | # Avec CI 3

# earth_GLM

#|n    |Methode |Erreur de Monte Carlo                  |Autres                                   |
#|:----|:-------|:--------------------------------------|:----------------------------------------|
#|1000 |TNDDR   |MCSE_bias         , 0.0451602906669299 |Bias               , -0.0739510311910888 |
#|-    |-       |MCSE_var          , 0.0105365036425814 |Var                , 0.020394518531216   |
#|-    |-       |MCSE_mse          , 0.0160031544285782 |Mse                , 0.0238238216923198  |
#|-    |-       |%Cov              , 0.9                |Précision_var      , 1.68006660404718    | # Avec CI 1
#|-    |-       |%Cov              , 0.9                |Précision_var      , 1.68006660404718    | # Avec CI 2
#|-    |-       |%Cov              , 0.9                |Précision_var      , 1.68006660404718    | # Avec CI 3

# Réseaux de neurones

#|n    |Methode |Erreur de Monte Carlo                  |Autres                                   |
#|:----|:-------|:--------------------------------------|:----------------------------------------|
#|1000 |TNDDR   |MCSE_bias        , 0.124471845358785   |Bias               , -0.0277326486617437 |
#|-    |-       |MCSE_var         , 0.0354189421907207  |Var                , 0.154932402870212   |
#|-    |-       |MCSE_mse         , 0.0326570246907765  |Mse                , 0.140208262384987   |
#|-    |-       |%Cov             , 1                   |Précision_var      , 1.68598960584829    | # Avec CI 1
#|-    |-       |%Cov             , 1                   |Précision_var      , 1.68598960584829    | # Avec CI 2
#|-    |-       |%Cov             , 1                   |Précision_var      , 1.68598960584829    | # Avec CI 3

