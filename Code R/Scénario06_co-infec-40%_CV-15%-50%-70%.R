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

################################################################################
# Scénario 06 (co-infection ~ 40% dans l'échantillon): Sous-scénario 01 : Couverture vaccinale = 15%, 50% et 70% 
#                                                                         I1_prev = 15% et 
#                                                                         I2 = 50%
# Pour CV = 15% : co_inf_para1 = 0.5, co_inf_para2 = 1
# Pour CV = 50% : co_inf_para1 = 0.5, co_inf_para2 = 0.5
# Pour CV = 70% : co_inf_para1 = 0.5, co_inf_para2 = 0.5

################################################################################

################################################################################
################################################################################
# Pour CV = 15% : co_inf_para1 = 0.5, co_inf_para2 = 1
# Proportion de co-infection ~40%

set.seed(1) # Pour avoir toujours les mêmes germes

nsim <- 10 # Nombre de réplications

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000, CV = 0.15,
                      co_inf_para1 = 0.5, co_inf_para2 = 1)
  
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

nsim <- 1000

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

resultats3 <- data.frame(matrix(ncol = 36, 
                                nrow = nsim))

colnames(resultats3) <- c("RRm_RF", "VE_RF", "var_log_RRm-RF", "IC_inf1-RF", "IC_sup1-RF", "IC_inf2-RF", "IC_sup2-RF", 
                          "IC_inf3-RF", "IC_sup3-RF",# Forêt aléatoire
                          
                          "RRm_Lasso", "VE_Lasso","var_log_RRm-Lasso", "IC_inf1-Lasso", "IC_sup1-Lasso", 
                          "IC_inf2-Lasso", "IC_sup2-Lasso", 
                          "IC_inf3-Lasso", "IC_sup3-Lasso",# Régression Lasso
                          
                          "RRm_Mars", "VE_Mars","var_log_RRm-Mars", "IC_inf1-Mars", "IC_sup1-Mars", 
                          "IC_inf2-Mars", "IC_sup2-Mars", 
                          "IC_inf3-Mars", "IC_sup3-Mars",# Régression avec splines
                          
                          "RRm_RN", "VE_RN","var_log_RRm-RN", "IC_inf1-RN", "IC_sup1-RN", 
                          "IC_inf2-RN", "IC_sup2-RN", 
                          "IC_inf3-RN", "IC_sup3-RN"# Réseaux de neurones
                          
                          
                          # Risque relatid marginal   
                          # Variance du log du risque relatif marginal
                          # Borne inférieure du premier intervalle de confiance
                          # Sa borne supérieure
                          # Borne inférieure du deuxième intervalle de confiance
                          # Sa borne supérieure
                          # Borne inférieure du troisième intervalle de confiance
                          # Sa borne supérieure
                          
)

methode <- list(RandomForest, Lasso, Mars, RN)

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, CV = 0.15,
                 co_inf_para1 = 0.5, co_inf_para2 = 1, popsize = 1*10**6)
  
  resultats[i,] <- RegLog(dat) # Régression logistique
  resultats2[i,] <- IPW(dat)   # IPW
  
  l <- list()
  
  for(j in methode) {
    k <- rep(NA, 9)
    tryCatch({
      
      k <- TNDDR(dat, j) # Liste des résultats pour la méthode j
      
    }, error = function(e){})
    
    l <- append(l, k) # Combiner les résultats des différentes méthodes
    
    resultats3[i,] <- l
    
    
  }
  
  
  # DT : Ajout d'une ligne pour suivre l'avancement
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}

# Combien de réplications avec des NA

sum(rowSums(is.na(resultats)) > 0)
sum(rowSums(is.na(resultats2)) > 0)

sum(rowSums(is.na(resultats3[,1:9])) > 0)
sum(rowSums(is.na(resultats3[,10:18])) > 0) #11
sum(rowSums(is.na(resultats3[,19:27])) > 0)
sum(rowSums(is.na(resultats3[,28:36])) > 0)
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
  list(name = "Var", value = MCSE_var[2]),
  list(name = "Mse", value = MCSE_mse[2]),
  list(name = "Précision_var", value = sd(resultats$coe_reg) - mean(resultats$err_reg)) # Voir si la variance est bien estimée
  
)

kable(Tab01)
# |n    |Methode |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:-------|:----------------------------------------|:--------------------------------------|
# |1000 |RegLog  |MCSE_bias          , 0.00527636962869153 |Bias             , 0.128320016455484   |
# |-    |-       |MCSE_var           , 0.00243943658910443 |Var               , 0.0278400764585784 |
# |-    |-       |MCSE_mse          , 0.0034131857855661   |Mse               , 0.0442782630052555 |
# |-    |-       |%Cov , 0.884                             |Précision_var     , 0.0114000582468003 |
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

Tab01$Methode <- c("IPW", "-", "-", "-")

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
# |n    |Methode |Erreur de Monte Carlo                  |Autres                                 |
# |:----|:-------|:--------------------------------------|:--------------------------------------|
# |1000 |IPW     |MCSE_bias        , 0.021047022627006   |Bias             , 0.156081830135965   |
# |-    |-       |MCSE_var         , 0.406220494715958   |Var              , 0.442977161461702   |
# |-    |-       |MCSE_mse          , 0.0034131857855661 |Mse               , 0.0442782630052555 |
# |-    |-       |%Cov, 0.9                              |Précision_var     , 0.0297091067074358 |

########################## TNDDR ###############################

Tab01$Methode <- c("TNDDR_RF", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_RF)
mean(resultats3$RRm_RF)
sd(resultats3$RRm_RF)

mean(sqrt(resultats3$`var_log_RRm-RF`))
sd(log(resultats3$RRm_RF))

help("calc_absolute") # calculer les différentes mesures de performance

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats3$vrai_param <- rep(mean(l_vraiRRm), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-RF`, `IC_sup2-RF`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RF)) - mean(sqrt(resultats3$`var_log_RRm-RF`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:--------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00652075825163294 |Bias              , 0.0113709629249469 |
# |-    |-        |MCSE_var           , 0.00214481368731423 |Var              , 0.042520288176239   |
# |-    |-        |MCSE_mse          , 0.0021695144284629   |Mse               , 0.0426070666859033 |
# |-    |-        |%Cov , 0.984                             |Précision_var    , 0.750627940473611   |

################################################################################

# Régression Lasso

Tab01$Methode <- c("TNDDR_Lasso", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_Lasso)
mean(resultats3$RRm_Lasso, na.rm = TRUE)
sd(resultats3$RRm_Lasso, na.rm = TRUE)

mean(sqrt(na.omit(resultats3$`var_log_RRm-Lasso`)))
sd(log(na.omit(resultats3$RRm_Lasso)))

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-Lasso`, `IC_sup2-Lasso`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(na.omit(resultats3$RRm_Lasso))) - mean(sqrt(na.omit(resultats3$`var_log_RRm-Lasso`))) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode     |Erreur de Monte Carlo                    |Autres                                   |
# |:----|:-----------|:----------------------------------------|:----------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias         , 0.0151030110730143   |Bias               , 0.00245904240709749 |
# |-    |-           |MCSE_var           , 0.00304987129022721 |Var              , 0.225591833093404     |
# |-    |-           |MCSE_mse           , 0.00310949791069763 |Mse              , 0.225369779039493     |
# |-    |-           |%Cov             , 0.981799797775531     |Précision_var    , -522.437389077505     |

################################################################################

# earth_GLM

Tab01$Methode <- c("TNDDR_Mars", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_Mars)
mean(resultats3$RRm_Mars)
sd(resultats3$RRm_Mars)

mean(sqrt(resultats3$`var_log_RRm-Mars`))
sd(log(resultats3$RRm_Mars))

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-Mars`, `IC_sup2-Mars`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Mars)) - mean(sqrt(resultats3$`var_log_RRm-Mars`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode    |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:----------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_Mars |MCSE_bias         , 0.0103706954552005   |Bias              , 0.0104678097370668 |
# |-    |-          |MCSE_var           , 0.00344984394539925 |Var              , 0.107551324224516   |
# |-    |-          |MCSE_mse           , 0.00352961746446772 |Mse              , 0.107553347940983   |
# |-    |-          |%Cov , 0.985                             |Précision_var    , -2.13377147344331   |

################################################################################

# Réseaux de neurones

Tab01$Methode <- c("TNDDR_RN", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_RN)
mean(resultats3$RRm_RN)
sd(resultats3$RRm_RN)

mean(sqrt(resultats3$`var_log_RRm-RN`))
sd(log(resultats3$RRm_RN))

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-RN`, `IC_sup2-RN`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RN)) - mean(sqrt(resultats3$`var_log_RRm-RN`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:--------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_RN |MCSE_bias         , 0.0114304770030273   |Bias             , -0.15563876371718   |
# |-    |-        |MCSE_var           , 0.00529711987744777 |Var              , 0.130655804516737   |
# |-    |-        |MCSE_mse           , 0.00290870000780567 |Mse              , 0.154748573483633   |
# |-    |-        |%Cov , 0.935                             |Précision_var     , -0.941195780347761 |
################################################################################
################################################################################
# Pour CV = 50% : co_inf_para1 = 1, co_inf_para2 = 2
# Proportion de co-infection ~40%

nsim <- 10 # Nombre de réplications

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000, CV = 0.5,
                      co_inf_para1 = 0.5, co_inf_para2 = 1)
  
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

nsim <- 1000

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

resultats3 <- data.frame(matrix(ncol = 36, 
                                nrow = nsim))

colnames(resultats3) <- c("RRm_RF", "VE_RF", "var_log_RRm-RF", "IC_inf1-RF", "IC_sup1-RF", "IC_inf2-RF", "IC_sup2-RF", 
                          "IC_inf3-RF", "IC_sup3-RF",# Forêt aléatoire
                          
                          "RRm_Lasso", "VE_Lasso","var_log_RRm-Lasso", "IC_inf1-Lasso", "IC_sup1-Lasso", 
                          "IC_inf2-Lasso", "IC_sup2-Lasso", 
                          "IC_inf3-Lasso", "IC_sup3-Lasso",# Régression Lasso
                          
                          "RRm_Mars", "VE_Mars","var_log_RRm-Mars", "IC_inf1-Mars", "IC_sup1-Mars", 
                          "IC_inf2-Mars", "IC_sup2-Mars", 
                          "IC_inf3-Mars", "IC_sup3-Mars",# Régression avec splines
                          
                          "RRm_RN", "VE_RN","var_log_RRm-RN", "IC_inf1-RN", "IC_sup1-RN", 
                          "IC_inf2-RN", "IC_sup2-RN", 
                          "IC_inf3-RN", "IC_sup3-RN"# Réseaux de neurones
                          
                          
                          # Risque relatid marginal   
                          # Variance du log du risque relatif marginal
                          # Borne inférieure du premier intervalle de confiance
                          # Sa borne supérieure
                          # Borne inférieure du deuxième intervalle de confiance
                          # Sa borne supérieure
                          # Borne inférieure du troisième intervalle de confiance
                          # Sa borne supérieure
                          
)

methode <- list(RandomForest, Lasso, Mars, RN)

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, CV = 0.5,
                 co_inf_para1 = 0.5, co_inf_para2 = 1, popsize = 1*10**6)
  
  resultats[i,] <- RegLog(dat) # Régression logistique
  resultats2[i,] <- IPW(dat)   # IPW
  
  l <- list()
  
  for(j in methode) {
    k <- rep(NA, 9)
    tryCatch({
      
      k <- TNDDR(dat, j) # Liste des résultats pour la méthode j
      
    }, error = function(e){})
    
    l <- append(l, k) # Combiner les résultats des différentes méthodes
    
    resultats3[i,] <- l
    
    
  }
  
  
  # DT : Ajout d'une ligne pour suivre l'avancement
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}

# Combien de réplications avec des NA

sum(rowSums(is.na(resultats)) > 0)
sum(rowSums(is.na(resultats2)) > 0)

sum(rowSums(is.na(resultats3[,1:9])) > 0)
sum(rowSums(is.na(resultats3[,10:18])) > 0)
sum(rowSums(is.na(resultats3[,19:27])) > 0)
sum(rowSums(is.na(resultats3[,28:36])) > 0)
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
  list(name = "Var", value = MCSE_var[2]),
  list(name = "Mse", value = MCSE_mse[2]),
  list(name = "Précision_var", value = sd(resultats$coe_reg) - mean(resultats$err_reg)) # Voir si la variance est bien estimée
  
)

kable(Tab01)
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-------|:------------------------------------------|:----------------------------------------|
# |1000 |RegLog  |MCSE_bias          , 0.00328996058953196   |Bias             , 0.121805714911302     |
# |-    |-       |MCSE_var            , 0.000560530534493368 |Var               , 0.0108238406806734   |
# |-    |-       |MCSE_mse          , 0.0011344005771718     |Mse               , 0.0256496490250462   |
# |-    |-       |%Cov , 0.734                               |Précision_var      , 0.00728250671968103 |
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

Tab01$Methode <- c("IPW", "-", "-", "-")

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
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                 |
# |:----|:-------|:------------------------------------------|:--------------------------------------|
# |1000 |IPW     |MCSE_bias          , 0.00342266705576123   |Bias              , 0.0958927950484479 |
# |-    |-       |MCSE_var            , 0.000670914584092675 |Var               , 0.0117146497745932 |
# |-    |-       |MCSE_mse          , 0.0011344005771718     |Mse               , 0.0256496490250462 |
# |-    |-       |%Cov, 0.85                                 |Précision_var     , -0.013289788098875 |

########################## TNDDR ###############################

Tab01$Methode <- c("TNDDR_RF", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_RF)
mean(resultats3$RRm_RF)
sd(resultats3$RRm_RF)

mean(sqrt(resultats3$`var_log_RRm-RF`))
sd(log(resultats3$RRm_RF))

help("calc_absolute") # calculer les différentes mesures de performance

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats3$vrai_param <- rep(mean(l_vraiRRm), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-RF`, `IC_sup2-RF`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RF)) - mean(sqrt(resultats3$`var_log_RRm-RF`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                      |Autres                                 |
# |:----|:--------|:------------------------------------------|:--------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00359444223570393   |Bias             , 0.105002962383612   |
# |-    |-        |MCSE_var            , 0.000725285512955618 |Var               , 0.0129200149858122 |
# |-    |-        |MCSE_mse           , 0.00119327203246633   |Mse               , 0.0239327170801607 |
# |-    |-        |%Cov , 0.948                               |Précision_var     , 0.0524540434100561 |

################################################################################

# Régression Lasso

Tab01$Methode <- c("TNDDR_Lasso", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_Lasso)
mean(resultats3$RRm_Lasso)
sd(resultats3$RRm_Lasso)

mean(sqrt(resultats3$`var_log_RRm-Lasso`))
sd(log(resultats3$RRm_Lasso))

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-Lasso`, `IC_sup2-Lasso`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Lasso)) - mean(sqrt(resultats3$`var_log_RRm-Lasso`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode     |Erreur de Monte Carlo                    |Autres                                   |
# |:----|:-----------|:----------------------------------------|:----------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias         , 0.0145110589726743   |Bias               , -0.0465776684032516 |
# |-    |-           |MCSE_var           , 0.00420385716458361 |Var              , 0.210570832508433     |
# |-    |-           |MCSE_mse           , 0.00307926210751053 |Mse              , 0.212529740869807     |
# |-    |-           |%Cov , 0.967                             |Précision_var    , -127.419832032829     |
################################################################################

# earth_GLM

Tab01$Methode <- c("TNDDR_Mars", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_Mars)
mean(resultats3$RRm_Mars)
sd(resultats3$RRm_Mars)

mean(sqrt(resultats3$`var_log_RRm-Mars`))
sd(log(resultats3$RRm_Mars))

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-Mars`, `IC_sup2-Mars`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Mars)) - mean(sqrt(resultats3$`var_log_RRm-Mars`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode    |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:----------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00781961179246053 |Bias              , 0.0522634612484528 |
# |-    |-          |MCSE_var           , 0.00282176030290267 |Var               , 0.0611463285847878 |
# |-    |-          |MCSE_mse           , 0.00285612419009826 |Mse               , 0.0638166516378716 |
# |-    |-          |%Cov , 0.919                             |Précision_var     , -0.783941870248265 |

################################################################################

# Réseaux de neurones

Tab01$Methode <- c("TNDDR_RN", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_RN)
mean(resultats3$RRm_RN)
sd(resultats3$RRm_RN)

mean(sqrt(resultats3$`var_log_RRm-RN`))
sd(log(resultats3$RRm_RN))

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-RN`, `IC_sup2-RN`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RN)) - mean(sqrt(resultats3$`var_log_RRm-RN`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                   |
# |:----|:--------|:----------------------------------------|:----------------------------------------|
# |1000 |TNDDR_RN |MCSE_bias         , 0.0112834142103624   |Bias               , -0.0816313216965155 |
# |-    |-        |MCSE_var           , 0.00417042245094675 |Var              , 0.127315436242608     |
# |-    |-        |MCSE_mse           , 0.00331068733505384 |Mse              , 0.133851793488285     |
# |-    |-        |%Cov , 0.792                             |Précision_var    , -5.64831507474623     |
################################################################################
################################################################################
# Pour CV = 70% : co_inf_para1 = 1, co_inf_para2 = 2
# Proportion de co-infection ~40%

nsim <- 10 # Nombre de réplications

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000, CV = 0.7,
                      co_inf_para1 = 0.5, co_inf_para2 = 1)
  
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

nsim <- 1000

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

resultats3 <- data.frame(matrix(ncol = 36, 
                                nrow = nsim))

colnames(resultats3) <- c("RRm_RF", "VE_RF", "var_log_RRm-RF", "IC_inf1-RF", "IC_sup1-RF", "IC_inf2-RF", "IC_sup2-RF", 
                          "IC_inf3-RF", "IC_sup3-RF",# Forêt aléatoire
                          
                          "RRm_Lasso", "VE_Lasso","var_log_RRm-Lasso", "IC_inf1-Lasso", "IC_sup1-Lasso", 
                          "IC_inf2-Lasso", "IC_sup2-Lasso", 
                          "IC_inf3-Lasso", "IC_sup3-Lasso",# Régression Lasso
                          
                          "RRm_Mars", "VE_Mars","var_log_RRm-Mars", "IC_inf1-Mars", "IC_sup1-Mars", 
                          "IC_inf2-Mars", "IC_sup2-Mars", 
                          "IC_inf3-Mars", "IC_sup3-Mars",# Régression avec splines
                          
                          "RRm_RN", "VE_RN","var_log_RRm-RN", "IC_inf1-RN", "IC_sup1-RN", 
                          "IC_inf2-RN", "IC_sup2-RN", 
                          "IC_inf3-RN", "IC_sup3-RN"# Réseaux de neurones
                          
                          
                          # Risque relatid marginal   
                          # Variance du log du risque relatif marginal
                          # Borne inférieure du premier intervalle de confiance
                          # Sa borne supérieure
                          # Borne inférieure du deuxième intervalle de confiance
                          # Sa borne supérieure
                          # Borne inférieure du troisième intervalle de confiance
                          # Sa borne supérieure
                          
)

methode <- list(RandomForest, Lasso, Mars, RN)

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, CV = 0.7,
                 co_inf_para1 = 0.5, co_inf_para2 = 1, popsize = 1*10**6)
  
  resultats[i,] <- RegLog(dat) # Régression logistique
  resultats2[i,] <- IPW(dat)   # IPW
  
  l <- list()
  
  for(j in methode) {
    k <- rep(NA, 9)
    tryCatch({
      
      k <- TNDDR(dat, j) # Liste des résultats pour la méthode j
      
    }, error = function(e){})
    
    l <- append(l, k) # Combiner les résultats des différentes méthodes
    
    resultats3[i,] <- l
    
    
  }
  
  
  # DT : Ajout d'une ligne pour suivre l'avancement
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}

# Combien de réplications avec des NA

sum(rowSums(is.na(resultats)) > 0)
sum(rowSums(is.na(resultats2)) > 0)

sum(rowSums(is.na(resultats3[,1:9])) > 0)
sum(rowSums(is.na(resultats3[,10:18])) > 0)
sum(rowSums(is.na(resultats3[,19:27])) > 0)
sum(rowSums(is.na(resultats3[,28:36])) > 0)
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
  list(name = "Var", value = MCSE_var[2]),
  list(name = "Mse", value = MCSE_mse[2]),
  list(name = "Précision_var", value = sd(resultats$coe_reg) - mean(resultats$err_reg)) # Voir si la variance est bien estimée
  
)

kable(Tab01)
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-------|:------------------------------------------|:----------------------------------------|
# |1000 |RegLog  |MCSE_bias          , 0.00347404157482489   |Bias             , 0.121871830570806     |
# |-    |-       |MCSE_var            , 0.000575355957048757 |Var               , 0.0120689648636118   |
# |-    |-       |MCSE_mse           , 0.00113658071168708   |Mse               , 0.0269096389854274   |
# |-    |-       |%Cov , 0.744                               |Précision_var      , 0.00866235592317038 |
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

Tab01$Methode <- c("IPW", "-", "-", "-")

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
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-------|:------------------------------------------|:----------------------------------------|
# |1000 |IPW     |MCSE_bias          , 0.00370033574810705   |Bias              , 0.0867333667928477   |
# |-    |-       |MCSE_var            , 0.000644165100534919 |Var               , 0.0136924846487189   |
# |-    |-       |MCSE_mse           , 0.00113658071168708   |Mse               , 0.0269096389854274   |
# |-    |-       |%Cov , 0.878                               |Précision_var      , -0.0200470567552394 |

########################## TNDDR ###############################

Tab01$Methode <- c("TNDDR_RF", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_RF)
mean(resultats3$RRm_RF)
sd(resultats3$RRm_RF)

mean(sqrt(resultats3$`var_log_RRm-RF`))
sd(log(resultats3$RRm_RF))

help("calc_absolute") # calculer les différentes mesures de performance

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats3$vrai_param <- rep(mean(l_vraiRRm), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-RF`, `IC_sup2-RF`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RF)) - mean(sqrt(resultats3$`var_log_RRm-RF`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:--------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00541802772001726 |Bias             , 0.138623152212134   |
# |-    |-        |MCSE_var           , 0.00164719080746032 |Var               , 0.0293550243748754 |
# |-    |-        |MCSE_mse           , 0.00251868211611588 |Mse               , 0.0485420476797291 |
# |-    |-        |%Cov , 0.948                             |Précision_var    , 0.285662619663488   |

################################################################################

# Régression Lasso

Tab01$Methode <- c("TNDDR_Lasso", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_Lasso)
mean(resultats3$RRm_Lasso)
sd(resultats3$RRm_Lasso)

mean(sqrt(resultats3$`var_log_RRm-Lasso`))
sd(log(resultats3$RRm_Lasso))

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-Lasso`, `IC_sup2-Lasso`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Lasso)) - mean(sqrt(resultats3$`var_log_RRm-Lasso`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode     |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:-----------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias         , 0.0138185854619369   |Bias              , -0.101842378306986 |
# |-    |-           |MCSE_var           , 0.00523201915156382 |Var              , 0.190953304168853   |
# |-    |-           |MCSE_mse           , 0.00285323645369409 |Mse              , 0.201134220883908   |
# |-    |-           |%Cov , 0.935                             |Précision_var    , -64.1572716012555   |

################################################################################

# earth_GLM

Tab01$Methode <- c("TNDDR_Mars", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_Mars)
mean(resultats3$RRm_Mars)
sd(resultats3$RRm_Mars)

mean(sqrt(resultats3$`var_log_RRm-Mars`))
sd(log(resultats3$RRm_Mars))

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-Mars`, `IC_sup2-Mars`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Mars)) - mean(sqrt(resultats3$`var_log_RRm-Mars`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode    |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:----------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00832403027596892 |Bias              , 0.0700714898226403 |
# |-    |-          |MCSE_var           , 0.00294866526935241 |Var               , 0.0692894800352473 |
# |-    |-          |MCSE_mse           , 0.00301311882659544 |Mse               , 0.0741302042411764 |
# |-    |-          |%Cov , 0.844                             |Précision_var     , 0.0554407116146658 |

################################################################################

# Réseaux de neurones

Tab01$Methode <- c("TNDDR_RN", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_RN)
mean(resultats3$RRm_RN)
sd(resultats3$RRm_RN)

mean(sqrt(resultats3$`var_log_RRm-RN`))
sd(log(resultats3$RRm_RN))

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- calc_coverage(resultats3, `IC_inf2-RN`, `IC_sup2-RN`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RN)) - mean(sqrt(resultats3$`var_log_RRm-RN`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:--------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_RN |MCSE_bias         , 0.0116492246073754   |Bias              , -0.109230698572112 |
# |-    |-        |MCSE_var           , 0.00451094831354888 |Var             , 0.13570443395308     |
# |-    |-        |MCSE_mse           , 0.00313796681842498 |Mse              , 0.147500075029679   |
# |-    |-        |%Cov , 0.639                             |Précision_var    , 0.135150762485929   |