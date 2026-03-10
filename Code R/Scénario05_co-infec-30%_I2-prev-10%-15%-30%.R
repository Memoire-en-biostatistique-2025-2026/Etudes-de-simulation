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
# Scénario 05 (co-infection ~ 40% dans l'échantillon): Sous-scénario 01 : I2_prev = 10%, 15% et 30% 
#                                                                         I1_prev = 15% et 
#                                                                         CV = 33%
# Pour I2_prev = 30% : co_inf_para1 = 3, co_inf_para2 = 5
# Pour I2_prev = 50% : co_inf_para1 = 3, co_inf_para2 = 5
# Pour I2_prev = 70% : co_inf_para1 = 1, co_inf_para2 = 0.5

# Pour I2_prev = 15% : co_inf_para1 = 4, co_inf_para2 = 5
################################################################################

################################################################################
################################################################################
# Pour I2_prev = 10% : co_inf_para1 = 3, co_inf_para2 = 5
# Proportion de co-infection ~30%

set.seed(1) # Pour avoir toujours les mêmes germes

nsim <- 10 # Nombre de réplications

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000, I2_prev = 0.1,
                      co_inf_para1 = 3, co_inf_para2 = 5)
  
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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, I2_prev = 0.1,
                 co_inf_para1 = 3, co_inf_para2 = 5, popsize = 1*10**6)
  
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
# |1000 |RegLog  |MCSE_bias          , 0.00276786122931905   |Bias             , 0.146580333616376     |
# |-    |-       |MCSE_var            , 0.000362705783498027 |Var                , 0.00766105578476757 |
# |-    |-       |MCSE_mse            , 0.000988982180119329 |Mse               , 0.0291391889320708   |
# |-    |-       |%Cov , 0.456                               |Précision_var    , 0.015901824854086     |
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
# |1000 |IPW     |MCSE_bias          , 0.00274493184269489   |Bias             , 0.140683166715318     |
# |-    |-       |MCSE_var            , 0.000354502330251623 |Var                , 0.00753465082104037 |
# |-    |-       |MCSE_mse            , 0.000988982180119329 |Mse               , 0.0291391889320708   |
# |-    |-       |%Cov, 0.41                                 |Précision_var     , 0.0310084121277522   |

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
# |n    |Methode  |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:--------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00286309600782421   |Bias             , 0.131984925642931     |
# |-    |-        |MCSE_var            , 0.000383675470712765 |Var                , 0.00819731875001892 |
# |-    |-        |MCSE_mse            , 0.000950867658322507 |Mse              , 0.025609142028239     |
# |-    |-        |%Cov , 0.702                               |Précision_var     , 0.0129262290262268   |

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
# |n    |Methode     |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-----------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias          , 0.00283254624754528   |Bias            , 0.15114624498264       |
# |-    |-           |MCSE_var            , 0.000373212598828243 |Var                , 0.00802331824448284 |
# |-    |-           |MCSE_mse           , 0.00103697514960632   |Mse               , 0.0308604822985907   |
# |-    |-           |%Cov , 0.326                               |Précision_var     , 0.0530995193175323   |

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
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00616703410096579 |Bias              , 0.0977401381554756 |
# |-    |-          |MCSE_var           , 0.00238660289474581 |Var               , 0.0380323096024749 |
# |-    |-          |MCSE_mse           , 0.00232758610459081 |Mse               , 0.0475474118995238 |
# |-    |-          |%Cov , 0.788                             |Précision_var   , 1.07046619786606     |

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
# |1000 |TNDDR_RN |MCSE_bias          , 0.00719901053249542 |Bias              , 0.0824470160310415 |
# |-    |-        |MCSE_var           , 0.00262810273038241 |Var             , 0.05182575264698     |
# |-    |-        |MCSE_mse           , 0.00257201955497939 |Mse               , 0.0585714373467559 |
# |-    |-        |%Cov , 0.826                             |Précision_var   , 1.31214499687336     |
################################################################################
################################################################################
# Pour I2_prev = 15% : co_inf_para1 = 3, co_inf_para2 = 5
# Proportion de co-infection ~30%

nsim <- 10 # Nombre de réplications

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000, I2_prev = 0.15,
                      co_inf_para1 = 3, co_inf_para2 = 5)
  
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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, I2_prev = 0.15,
                 co_inf_para1 = 3, co_inf_para2 = 5, popsize = 1*10**6)
  
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
# |1000 |RegLog  |MCSE_bias          , 0.00268531236999125   |Bias             , 0.147952039385118     |
# |-    |-       |MCSE_var            , 0.000369012943829343 |Var                , 0.00721090252442803 |
# |-    |-       |MCSE_mse            , 0.000984230179910968 |Mse               , 0.0290934975801191   |
# |-    |-       |%Cov , 0.459                               |Précision_var     , 0.0102191886164966   |
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
# |1000 |IPW     |MCSE_bias          , 0.00267982963504437   |Bias             , 0.138484557808066     |
# |-    |-       |MCSE_var            , 0.000361540853137267 |Var                , 0.00718148687286205 |
# |-    |-       |MCSE_mse            , 0.000984230179910968 |Mse               , 0.0290934975801191   |
# |-    |-       |%Cov , 0.453                               |Précision_var     , 0.0233409419464073   |

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
# |n    |Methode  |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:--------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00281962807532733   |Bias             , 0.127183466883913     |
# |-    |-        |MCSE_var            , 0.000387088061701697 |Var                , 0.00795030248317411 |
# |-    |-        |MCSE_mse            , 0.000907077817051087 |Mse               , 0.0241179864293024   |
# |-    |-        |%Cov , 0.723                               |Précision_var      , 0.00895027310176214 |

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
# |n    |Methode     |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-----------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias          , 0.00296261643027887   |Bias             , 0.150573282172128     |
# |-    |-           |MCSE_var            , 0.000679391092881278 |Var                , 0.00877709611295831 |
# |-    |-           |MCSE_mse           , 0.00112233401570415   |Mse               , 0.0314406323209326   |
# |-    |-           |%Cov , 0.273                               |Précision_var    , 0.278009310276432     |
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
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00624340767862687 |Bias              , 0.0977461928138397 |
# |-    |-          |MCSE_var           , 0.00241700796786513 |Var               , 0.0389801394415369 |
# |-    |-          |MCSE_mse           , 0.00241908157946957 |Mse               , 0.0484954775116957 |
# |-    |-          |%Cov , 0.828                             |Précision_var    , 0.956626735905759   |

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
################################################################################
################################################################################
# Pour I2_prev = 30% : co_inf_para1 = 1, co_inf_para2 = 0.5
# Proportion de co-infection ~30%

nsim <- 10 # Nombre de réplications

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000, I2_prev = 0.3,
                      co_inf_para1 = 1, co_inf_para2 = 0.5)
  
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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, I2_prev = 0.3,
                 co_inf_para1 = 1, co_inf_para2 = 0.5, popsize = 1*10**6)
  
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
# |1000 |RegLog  |MCSE_bias          , 0.00284729194753696   |Bias             , 0.108008237386338     |
# |-    |-       |MCSE_var            , 0.000406168959768163 |Var                , 0.00810707143450883 |
# |-    |-       |MCSE_mse            , 0.000856335764173157 |Mse               , 0.0197647437063779   |
# |-    |-       |%Cov , 0.685                               |Précision_var     , 0.0159674041820914   |
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
# |1000 |IPW     |MCSE_bias          , 0.00282011440138918   |Bias              , 0.0938620062555732   |
# |-    |-       |MCSE_var            , 0.000397793242065879 |Var                , 0.00795304523692265 |
# |-    |-       |MCSE_mse            , 0.000856335764173157 |Mse               , 0.0197647437063779   |
# |-    |-       |%Cov , 0.739                               |Précision_var     , 0.0174025677651279   |

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
# |n    |Methode  |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:--------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00294368955265658   |Bias              , 0.0784413662397164   |
# |-    |-        |MCSE_var            , 0.000417803287444953 |Var                , 0.00866530818241947 |
# |-    |-        |MCSE_mse            , 0.000716423854539482 |Mse               , 0.0148096908117904   |
# |-    |-        |%Cov, 0.91                                 |Précision_var      , 0.00893187432132136 |

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
# |1000 |TNDDR_Lasso |MCSE_bias         , 0.0072245514447218   |Bias             , 0.125169476540054   |
# |-    |-           |MCSE_var           , 0.00294967089254731 |Var               , 0.0521941435774318 |
# |-    |-           |MCSE_mse           , 0.00334703897225405 |Mse               , 0.0678093472911656 |
# |-    |-           |%Cov , 0.477                             |Précision_var   , 1.14083839011641     |
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
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00662483674700934 |Bias              , 0.0454222949558141 |
# |-    |-          |MCSE_var           , 0.00246743211200248 |Var               , 0.0438884619245252 |
# |-    |-          |MCSE_mse           , 0.00250437350126885 |Mse               , 0.0459077583416537 |
# |-    |-          |%Cov , 0.945                             |Précision_var   , 1.08513872967487     |

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
# |1000 |TNDDR_RN |MCSE_bias          , 0.00914879323558243 |Bias               , -0.0231870280784262 |
# |-    |-        |MCSE_var          , 0.0030704244029703   |Var               , 0.0837004176674388   |
# |-    |-        |MCSE_mse           , 0.00301527942150028 |Mse               , 0.0841543555208811   |
# |-    |-        |%Cov , 0.901                             |Précision_var    , 0.884704420487132     |