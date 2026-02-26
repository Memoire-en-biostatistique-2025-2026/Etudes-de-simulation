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
# Scénario 03 (co-infection ~ 40% dans l'échantillon): Sous-scénario 01 : I2_prev = 30%, 50% et 70% 
#                                                                         I1_prev = 15% et 
#                                                                         CV = 33%
# Pour I2_prev = 30% : co_inf_para1 = 1, co_inf_para2 = 2
# Pour I2_prev = 50% : co_inf_para1 = 1, co_inf_para2 = 2
# Pour I2_prev = 70% : co_inf_para1 = 1, co_inf_para2 = 2

# Pour I2_prev = 15% : co_inf_para1 = 4, co_inf_para2 = 5
################################################################################

# Listes des valeurs à tester (à discuter)

I2_prev_liste <- c(0.3, 0.5, 0.7)

# Calcul des vraies valeurs 

set.seed(1) # Pour avoir toujours les mêmes germes

nsim1 <- 10 # Nombre de réplications

seeds_list <- sample(1:1000000, size = nsim1*3)

l_vraiRRc <- rep(NA, nsim1*3)
l_vraiRRm <- rep(NA, nsim1*3)

d <- 1
f <- nsim1

for(j in I2_prev_liste) {
  
  for (i in d:f) {
    
    dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000,
                        co_inf_para1 = 1, co_inf_para2 = 2, I2_prev = j)
    
    summary(dat)
    
    dat0 <- data.frame(C = dat$C, V = 0, Y = dat$I2_0*dat$W2_0*dat$H_0)
    
    dat1 <- data.frame(C = dat$C, V = 1, Y = dat$I2_1*dat$W2_1*dat$H_1) 
    
    dat_complet <- rbind(dat0, dat1)
    
    vraiRRc <- glm(Y ~ V + C, family = binomial(link = "logit"), data = dat_complet)
    
    l_vraiRRc[i] <- exp(coef(vraiRRc)[2])
    
    l_vraiRRm[i] <- mean(dat1$Y)/mean(dat0$Y)
    
    print(data.frame(Sys.time(), i))
  }
  
  d <- d + nsim1
  f <- f + nsim1
  
}

Tab <- data.frame( 
  
  `Prévalence de I2` = c(0.3, 0.5, 0.7),
  Mean_l_vraiRRc = c(mean(l_vraiRRc[1:nsim1]), mean(l_vraiRRc[(nsim1 + 1):(2*nsim1)]), mean(l_vraiRRc[(2*nsim1 + 1):(3*nsim1)])) ,
  SD_l_vraiRRc = c(sd(l_vraiRRc[1:nsim1]), sd(l_vraiRRc[(nsim1 + 1):(2*nsim1)]),sd(l_vraiRRc[(2*nsim1 + 1):(3*nsim1)])),
  Mean_l_vraiRRm =  c(mean(l_vraiRRm[1:nsim1]), mean(l_vraiRRm[(nsim1 + 1):(2*nsim1)]), mean(l_vraiRRm[(2*nsim1 + 1):(3*nsim1)])),
  SD_l_vraiRRm = c(sd(l_vraiRRm[1:nsim1]), sd(l_vraiRRm[(nsim1 + 1):(2*nsim1)]), sd(l_vraiRRm[(2*nsim1 + 1):(3*nsim1)]))
  
)

kable(Tab)

################################################################################
# Initialiser des objets pour contenir les resultats

Tab01 <- data.frame(n = c("5000", "-", "-", "-"))

################################################################################
########################## Analyse des résultats ###############################

nsim <- 1000

# Initialiser des objets pour contenir les resultats

resultats <- data.frame(matrix(ncol = 4, 
                               nrow = nsim * 3))

colnames(resultats) <- c("coe_reg", 
                         "err_reg",
                         "RRc", # Risque relatif conditionnel
                         "est_VE")

resultats2 <- data.frame(matrix(ncol = 5, 
                                nrow = nsim * 3))

colnames(resultats2) <- c("RRm",# Risque relatif marginal
                          "VE",
                          "var_log_RRm",# Variance du log du risque relatif marginal
                          "IC_inf", # Borne inférieure de l'intervalle de confiance
                          "IC_sup") # Sa borne supérieure

resultats3 <- data.frame(matrix(ncol = 36, 
                                nrow = nsim * 3))

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

############################### Estimation #####################################

d <- 1
f <- nsim

for(j in I2_prev_liste) {
  
  for (i in d:f) {
    
    
    dat <- datagen(seed = seeds_list[i], ssize = 5000, co_inf_para1 = 1, co_inf_para2 = 2, I2_prev = j, popsize = 1*10**6)
    
    resultats[i,] <- RegLog(dat) # Régression logistique
    resultats2[i,] <- IPW(dat)   # IPW
    
    l <- list()
    
    for(k in methode) {
      
      m <- TNDDR(dat, k) # Liste des résultats pour la méthode k
      l <- append(l, m) # Combiner les résultats des différentes méthodes
      
    }
    
    resultats3[i,] <- l
    
    
    # DT : Ajout d'une ligne pour suivre l'avancement
    
    if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
    
  }
  
  d <- d + nsim
  f <- f + nsim
  
}


########################## Affichage des résultats #############################

resultats$`Prévalence-I2` <- c("0.3", rep("-", nsim - 1,), "0.5", rep("-", nsim - 1,), "0.7", rep("-", nsim - 1,))
resultats$`Prévalence-I1` <- c("0.15", rep("-", 3*nsim - 1,))
resultats$`Couverture vaccinale` <- c("0.33", rep("-", 3*nsim - 1,))

resultats2$`Prévalence-I2` <- c("0.3", rep("-", nsim - 1,), "0.5", rep("-", nsim - 1,), "0.7", rep("-", nsim - 1,))
resultats2$`Prévalence-I1` <- c("0.15", rep("-", 3*nsim - 1,))
resultats2$`Couverture vaccinale` <- c("0.33", rep("-", 3*nsim - 1,))


resultats3$`Prévalence-I2` <- c("0.3", rep("-", nsim - 1,), "0.5", rep("-", nsim - 1,), "0.7", rep("-", nsim - 1,))
resultats3$`Prévalence-I1` <- c("0.15", rep("-", 3*nsim - 1,))
resultats3$`Couverture vaccinale` <- c("0.33", rep("-", 3*nsim - 1,))

################################################################################

########################## Régression logistique ###############################

##    - statistiques descriptives/valeur de couverture vaccinale

Tab <- data.frame( 
  
  n = c("5000", "-", "-"), Methode = c("RegLog", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats$RRc[1:nsim]), mean(resultats$RRc[(nsim + 1):(2*nsim)]), mean(resultats$RRc[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats$RRc[1:nsim]), sd(resultats$RRc[(nsim + 1):(2*nsim)]),sd(resultats$RRc[(2*nsim + 1):(3*nsim)])),
  Mean_err_reg =  c(mean(resultats$err_reg[1:nsim]), mean(resultats$err_reg[(nsim + 1):(2*nsim)]), mean(resultats$err_reg[(2*nsim + 1):(3*nsim)])),
  SD_coe_reg = c(sd(resultats$coe_reg[1:nsim]), sd(resultats$coe_reg[(nsim + 1):(2*nsim)]), sd(resultats$coe_reg[(2*nsim + 1):(3*nsim)]))
  
)

kable(Tab)

help("calc_absolute") # calculer les différentes mesures de performance

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats$vrai_param <- rep(c(mean(l_vraiRRc[1:nsim1]), mean(l_vraiRRc[(nsim1 + 1):(2*nsim1)]), mean(l_vraiRRc[(2*nsim1 + 1):(3*nsim1)])), each = nsim)

### Biais et MCSE_biais

Biais <- c(
  
  calc_absolute(resultats[1:nsim,], RRc, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats[(nsim + 1):(2*nsim),], RRc, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats[(2*nsim + 1):(3*nsim),], RRc, vrai_param, criteria = "bias")[2]
  
)

MCSE_biais <- c(
  
  calc_absolute(resultats[1:nsim,], RRc, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats[(nsim + 1):(2*nsim),], RRc, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats[(2*nsim + 1):(3*nsim),], RRc, vrai_param, criteria = "bias")[3]
  
)

### Var et MCSE_var

Var <- c(
  
  calc_absolute(resultats[1:nsim,], RRc, vrai_param, criteria = "var")[2],
  calc_absolute(resultats[(nsim + 1):(2*nsim),], RRc, vrai_param, criteria = "var")[2],
  calc_absolute(resultats[(2*nsim + 1):(3*nsim),], RRc, vrai_param, criteria = "var")[2]
  
)

MCSE_var <- c(
  
  calc_absolute(resultats[1:nsim,], RRc, vrai_param, criteria = "var")[3],
  calc_absolute(resultats[(nsim + 1):(2*nsim),], RRc, vrai_param, criteria = "var")[3],
  calc_absolute(resultats[(2*nsim + 1):(3*nsim),], RRc, vrai_param, criteria = "var")[3]
  
)

### MSE et MCSE_MSE 

Mse <- c(
  
  calc_absolute(resultats[1:nsim,], RRc, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats[(nsim + 1):(2*nsim),], RRc, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats[(2*nsim + 1):(3*nsim),], RRc, vrai_param, criteria = "mse")[2]
  
)

MCSE_mse <- c(
  
  calc_absolute(resultats[1:nsim,], RRc, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats[(nsim + 1):(2*nsim),], RRc, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats[(2*nsim + 1):(3*nsim),], RRc, vrai_param, criteria = "mse")[3]
  
)

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

resultats$lim_inf <- c(
  
  exp(resultats[1:nsim, 1] - 1.96*resultats[1:nsim, 2]), 
  exp(resultats[(nsim + 1):(2*nsim), 1] - 1.96*resultats[(nsim + 1):(2*nsim), 2]),
  exp(resultats[(2*nsim + 1):(3*nsim), 1] - 1.96*resultats[(2*nsim + 1):(3*nsim), 2])
  
)

resultats$lim_sup <- c(
  
  exp(resultats[1:nsim, 1] + 1.96*resultats[1:nsim, 2]), 
  exp(resultats[(nsim + 1):(2*nsim), 1] + 1.96*resultats[(nsim + 1):(2*nsim), 2]),
  exp(resultats[(2*nsim + 1):(3*nsim), 1] + 1.96*resultats[(2*nsim + 1):(3*nsim), 2])
  
)

coverage <- c(
  
  calc_coverage(resultats[1:nsim,], lim_inf, lim_sup, vrai_param)[2],
  calc_coverage(resultats[(nsim + 1):(2*nsim), ], lim_inf, lim_sup, vrai_param)[2],
  calc_coverage(resultats[(2*nsim + 1):(3*nsim), ], lim_inf, lim_sup, vrai_param)[2]
  
)

Precision_var <- c(
  
  sd(resultats$coe_reg[1:nsim]) - mean(resultats$err_reg[1:nsim]),
  sd(resultats$coe_reg[(nsim + 1):(2*nsim)]) - mean(resultats$err_reg[(nsim + 1):(2*nsim)]),
  sd(resultats$coe_reg[(2*nsim + 1):(3*nsim)]) - mean(resultats$err_reg[(2*nsim + 1):(3*nsim)])
  
)

Tab01$`Erreur de Monte Carlo_(30%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(30%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)

################################## IPW #########################################

##    - statistiques descriptives/valeur de couverture vaccinale

Tab <- data.frame( 
  
  n = c("5000", "-", "-"), Methode = c("IPW", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats2$RRm[1:nsim]), mean(resultats2$RRm[(nsim + 1):(2*nsim)]), mean(resultats2$RRm[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats2$RRm[1:nsim]), sd(resultats2$RRm[(nsim + 1):(2*nsim)]),sd(resultats2$RRm[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats2$var_log_RRm[1:nsim])), mean(sqrt(resultats2$var_log_RRm[(nsim + 1):(2*nsim)])), mean(sqrt(resultats2$var_log_RRm[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats2$RRm[1:nsim])), sd(log(resultats2$RRm[(nsim + 1):(2*nsim)])), sd(log(resultats2$RRm[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)

help("calc_absolute") # calculer les différentes mesures de performance

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats2$vrai_param <- rep(c(mean(l_vraiRRm[1:nsim1]), mean(l_vraiRRm[(nsim1 + 1):(2*nsim1)]), mean(l_vraiRRm[(2*nsim1 + 1):(3*nsim1)])), each = nsim)

### Biais et MCSE_biais

Biais <- c(
  
  calc_absolute(resultats2[1:nsim,], RRm, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats2[(nsim + 1):(2*nsim),], RRm, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats2[(2*nsim + 1):(3*nsim),], RRm, vrai_param, criteria = "bias")[2]
  
)

MCSE_biais <- c(
  
  calc_absolute(resultats2[1:nsim,], RRm, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats2[(nsim + 1):(2*nsim),], RRm, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats2[(2*nsim + 1):(3*nsim),], RRm, vrai_param, criteria = "bias")[3]
  
)

### Var et MCSE_var

Var <- c(
  
  calc_absolute(resultats2[1:nsim,], RRm, vrai_param, criteria = "var")[2],
  calc_absolute(resultats2[(nsim + 1):(2*nsim),], RRm, vrai_param, criteria = "var")[2],
  calc_absolute(resultats2[(2*nsim + 1):(3*nsim),], RRm, vrai_param, criteria = "var")[2]
  
)

MCSE_var <- c(
  
  calc_absolute(resultats2[1:nsim,], RRm, vrai_param, criteria = "var")[3],
  calc_absolute(resultats2[(nsim + 1):(2*nsim),], RRm, vrai_param, criteria = "var")[3],
  calc_absolute(resultats2[(2*nsim + 1):(3*nsim),], RRm, vrai_param, criteria = "var")[3]
  
)

### MSE et MCSE_MSE 

Mse <- c(
  
  calc_absolute(resultats2[1:nsim,], RRm, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats2[(nsim + 1):(2*nsim),], RRm, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats2[(2*nsim + 1):(3*nsim),], RRm, vrai_param, criteria = "mse")[2]
  
)

MCSE_mse <- c(
  
  calc_absolute(resultats2[1:nsim,], RRm, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats2[(nsim + 1):(2*nsim),], RRm, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats2[(2*nsim + 1):(3*nsim),], RRm, vrai_param, criteria = "mse")[3]
  
)

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- c(
  
  calc_coverage(resultats2[1:nsim,], IC_inf, IC_sup, vrai_param)[2],
  calc_coverage(resultats2[(nsim + 1):(2*nsim), ], IC_inf, IC_sup, vrai_param)[2],
  calc_coverage(resultats2[(2*nsim + 1):(3*nsim), ], IC_inf, IC_sup, vrai_param)[2]
  
)

Precision_var <- c(
  
  sd(log(resultats2$RRm[1:nsim])) - mean(sqrt(resultats2$var_log_RRm[1:nsim])),
  sd(log(resultats2$RRm[(nsim + 1):(2*nsim)])) - mean(sqrt(resultats2$var_log_RRm[(nsim + 1):(2*nsim)])),
  sd(log(resultats2$RRm[(2*nsim + 1):(3*nsim)])) - mean(sqrt(resultats2$var_log_RRm[(2*nsim + 1):(3*nsim)]))
  
)

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

Tab01$`Erreur de Monte Carlo_(30%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(30%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)


################################## TNDDR #######################################

# Avec forêt aléatoire

Tab <- data.frame( 
  
  n = c("5000", "-", "-"), Methode = c("TNDDR_FR", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_RF[1:nsim]), mean(resultats3$RRm_RF[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_RF[1:nsim]), sd(resultats3$RRm_RF[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-RF`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-RF`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-RF`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_RF[1:nsim])), sd(log(resultats3$RRm_RF[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)

help("calc_absolute") # calculer les différentes mesures de performance

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats3$vrai_param <- rep(c(mean(l_vraiRRm[1:nsim1]), mean(l_vraiRRm[(nsim1 + 1):(2*nsim1)]), mean(l_vraiRRm[(2*nsim1 + 1):(3*nsim1)])), each = nsim)

### Biais et MCSE_biais

Biais <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RF, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RF, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RF, vrai_param, criteria = "bias")[2]
  
)

MCSE_biais <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RF, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RF, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RF, vrai_param, criteria = "bias")[3]
  
)

### Var et MCSE_var

Var <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RF, vrai_param, criteria = "var")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RF, vrai_param, criteria = "var")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RF, vrai_param, criteria = "var")[2]
  
)

MCSE_var <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RF, vrai_param, criteria = "var")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RF, vrai_param, criteria = "var")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RF, vrai_param, criteria = "var")[3]
  
)

### MSE et MCSE_MSE 

Mse <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RF, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RF, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RF, vrai_param, criteria = "mse")[2]
  
)

MCSE_mse <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RF, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RF, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RF, vrai_param, criteria = "mse")[3]
  
)

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- c(
  
  calc_coverage(resultats3[1:nsim,], `IC_inf2-RF`, `IC_sup2-RF`, vrai_param)[2],
  calc_coverage(resultats3[(nsim + 1):(2*nsim), ], `IC_inf2-RF`, `IC_sup2-RF`, vrai_param)[2],
  calc_coverage(resultats3[(2*nsim + 1):(3*nsim), ], `IC_inf2-RF`, `IC_sup2-RF`, vrai_param)[2]
  
)

Precision_var <- c(
  
  sd(log(resultats3$RRm_RF[1:nsim])) - mean(sqrt(resultats3$`var_log_RRm-RF`[1:nsim])),
  sd(log(resultats3$RRm_RF[(nsim + 1):(2*nsim)])) - mean(sqrt(resultats3$`var_log_RRm-RF`[(nsim + 1):(2*nsim)])),
  sd(log(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])) - mean(sqrt(resultats3$`var_log_RRm-RF`[(2*nsim + 1):(3*nsim)]))
  
)

Tab01 <- data.frame(n = c("5000", "-", "-", "-"))

Tab01$`Erreur de Monte Carlo_(30%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(30%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)

################################################################################
# Régression Lasso

Tab <- data.frame( 
  
  n = c("5000", "-", "-"), Methode = c("TNDDR_Lasso", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_Lasso[1:nsim]), mean(resultats3$RRm_Lasso[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_Lasso[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_Lasso[1:nsim]), sd(resultats3$RRm_Lasso[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_Lasso[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-Lasso`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-Lasso`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-Lasso`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_Lasso[1:nsim])), sd(log(resultats3$RRm_Lasso[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_Lasso[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)


help("calc_absolute") # calculer les différentes mesures de performance

### Biais et MCSE_biais

Biais <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Lasso, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Lasso, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Lasso, vrai_param, criteria = "bias")[2]
  
)

MCSE_biais <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Lasso, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Lasso, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Lasso, vrai_param, criteria = "bias")[3]
  
)

### Var et MCSE_var

Var <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Lasso, vrai_param, criteria = "var")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Lasso, vrai_param, criteria = "var")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Lasso, vrai_param, criteria = "var")[2]
  
)

MCSE_var <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Lasso, vrai_param, criteria = "var")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Lasso, vrai_param, criteria = "var")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Lasso, vrai_param, criteria = "var")[3]
  
)

### MSE et MCSE_MSE 

Mse <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Lasso, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Lasso, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Lasso, vrai_param, criteria = "mse")[2]
  
)

MCSE_mse <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Lasso, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Lasso, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Lasso, vrai_param, criteria = "mse")[3]
  
)

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- c(
  
  calc_coverage(resultats3[1:nsim,], `IC_inf2-Lasso`, `IC_sup2-Lasso`, vrai_param)[2],
  calc_coverage(resultats3[(nsim + 1):(2*nsim), ], `IC_inf2-Lasso`, `IC_sup2-Lasso`, vrai_param)[2],
  calc_coverage(resultats3[(2*nsim + 1):(3*nsim), ], `IC_inf2-Lasso`, `IC_sup2-Lasso`, vrai_param)[2]
  
)

Precision_var <- c(
  
  sd(log(resultats3$RRm_RF[1:nsim])) - mean(sqrt(resultats3$`var_log_RRm-Lasso`[1:nsim])),
  sd(log(resultats3$RRm_RF[(nsim + 1):(2*nsim)])) - mean(sqrt(resultats3$`var_log_RRm-Lasso`[(nsim + 1):(2*nsim)])),
  sd(log(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])) - mean(sqrt(resultats3$`var_log_RRm-Lasso`[(2*nsim + 1):(3*nsim)]))
  
)

Tab01 <- data.frame(n = c("5000", "-", "-", "-"))

Tab01$`Erreur de Monte Carlo_(30%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(30%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)

################################################################################
# Régression avec splines: MARS

Tab <- data.frame( 
  
  n = c("5000", "-", "-"), Methode = c("TNDDR_Mars", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_Mars[1:nsim]), mean(resultats3$RRm_Mars[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_Mars[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_Mars[1:nsim]), sd(resultats3$RRm_Mars[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_Mars[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-Mars`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-Mars`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-Mars`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_Mars[1:nsim])), sd(log(resultats3$RRm_Mars[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_Mars[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)

help("calc_absolute") # calculer les différentes mesures de performance

### Biais et MCSE_biais

Biais <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Mars, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Mars, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Mars, vrai_param, criteria = "bias")[2]
  
)

MCSE_biais <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Mars, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Mars, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Mars, vrai_param, criteria = "bias")[3]
  
)

### Var et MCSE_var

Var <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Mars, vrai_param, criteria = "var")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Mars, vrai_param, criteria = "var")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Mars, vrai_param, criteria = "var")[2]
  
)

MCSE_var <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Mars, vrai_param, criteria = "var")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Mars, vrai_param, criteria = "var")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Mars, vrai_param, criteria = "var")[3]
  
)

### MSE et MCSE_MSE 

Mse <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Mars, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Mars, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Mars, vrai_param, criteria = "mse")[2]
  
)

MCSE_mse <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_Mars, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_Mars, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_Mars, vrai_param, criteria = "mse")[3]
  
)

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- c(
  
  calc_coverage(resultats3[1:nsim,], `IC_inf2-Mars`, `IC_sup2-Mars`, vrai_param)[2],
  calc_coverage(resultats3[(nsim + 1):(2*nsim), ], `IC_inf2-Mars`, `IC_sup2-Mars`, vrai_param)[2],
  calc_coverage(resultats3[(2*nsim + 1):(3*nsim), ], `IC_inf2-Mars`, `IC_sup2-Mars`, vrai_param)[2]
  
)

Precision_var <- c(
  
  sd(log(resultats3$RRm_Mars[1:nsim])) - mean(sqrt(resultats3$`var_log_RRm-Mars`[1:nsim])),
  sd(log(resultats3$RRm_Mars[(nsim + 1):(2*nsim)])) - mean(sqrt(resultats3$`var_log_RRm-Mars`[(nsim + 1):(2*nsim)])),
  sd(log(resultats3$RRm_Mars[(2*nsim + 1):(3*nsim)])) - mean(sqrt(resultats3$`var_log_RRm-Mars`[(2*nsim + 1):(3*nsim)]))
  
)

Tab01 <- data.frame(n = c("5000", "-", "-", "-"))

Tab01$`Erreur de Monte Carlo_(30%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(30%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)

################################################################################
# Réseaux de neurones

Tab <- data.frame( 
  
  n = c("5000", "-", "-"), Methode = c("TNDDR_RN", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_RN[1:nsim]), mean(resultats3$RRm_RN[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_RN[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_RN[1:nsim]), sd(resultats3$RRm_RN[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_RN[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-RN`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-RN`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-RN`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_RN[1:nsim])), sd(log(resultats3$RRm_RN[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_RN[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)

help("calc_absolute") # calculer les différentes mesures de performance

### Biais et MCSE_biais

Biais <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RN, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RN, vrai_param, criteria = "bias")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RN, vrai_param, criteria = "bias")[2]
  
)

MCSE_biais <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RN, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RN, vrai_param, criteria = "bias")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RN, vrai_param, criteria = "bias")[3]
  
)

### Var et MCSE_var

Var <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RN, vrai_param, criteria = "var")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RN, vrai_param, criteria = "var")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RN, vrai_param, criteria = "var")[2]
  
)

MCSE_var <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RN, vrai_param, criteria = "var")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RN, vrai_param, criteria = "var")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RN, vrai_param, criteria = "var")[3]
  
)

### MSE et MCSE_MSE 

Mse <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RN, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RN, vrai_param, criteria = "mse")[2],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RN, vrai_param, criteria = "mse")[2]
  
)

MCSE_mse <- c(
  
  calc_absolute(resultats3[1:nsim,], RRm_RN, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats3[(nsim + 1):(2*nsim),], RRm_RN, vrai_param, criteria = "mse")[3],
  calc_absolute(resultats3[(2*nsim + 1):(3*nsim),], RRm_RN, vrai_param, criteria = "mse")[3]
  
)

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

coverage <- c(
  
  calc_coverage(resultats3[1:nsim,], `IC_inf2-RN`, `IC_sup2-RN`, vrai_param)[2],
  calc_coverage(resultats3[(nsim + 1):(2*nsim), ], `IC_inf2-RN`, `IC_sup2-RN`, vrai_param)[2],
  calc_coverage(resultats3[(2*nsim + 1):(3*nsim), ], `IC_inf2-RN`, `IC_sup2-RN`, vrai_param)[2]
  
)

Precision_var <- c(
  
  sd(log(resultats3$RRm_RF[1:nsim])) - mean(sqrt(resultats3$`var_log_RRm-RN`[1:nsim])),
  sd(log(resultats3$RRm_RF[(nsim + 1):(2*nsim)])) - mean(sqrt(resultats3$`var_log_RRm-RN`[(nsim + 1):(2*nsim)])),
  sd(log(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])) - mean(sqrt(resultats3$`var_log_RRm-RN`[(2*nsim + 1):(3*nsim)]))
  
)

Tab01 <- data.frame(n = c("5000", "-", "-", "-"))

Tab01$`Erreur de Monte Carlo_(30%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(30%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)

################################################################################
################################################################################
# Pour I2_prev = 15% : co_inf_para1 = 4, co_inf_para2 = 5
# Proportion de co-infection ~36%

nsim <- 10 # Nombre de réplications

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000, I2_prev = 0.15,
                      co_inf_para1 = 4, co_inf_para2 = 5)
  
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

Tab01 <- data.frame(n = c("5000", "-", "-", "-"))

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
  
  dat <- datagen(seed = seeds_list[i], ssize = 5000, I2_prev = 0.15,
                 co_inf_para1 = 4, co_inf_para2 = 5, popsize = 1*10**6)
  
  resultats[i,] <- RegLog(dat) # Régression logistique
  resultats2[i,] <- IPW(dat)   # IPW
  
  l <- list()
  
  for(j in methode) {
    
    k <- TNDDR(dat, j) # Liste des résultats pour la méthode j
    l <- append(l, k) # Combiner les résultats des différentes méthodes
    
  }
  
  resultats3[i,] <- l
  
  
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
  list(name = "Var", value = MCSE_var[2]),
  list(name = "Mse", value = MCSE_mse[2]),
  list(name = "Précision_var", value = sd(resultats$coe_reg) - mean(resultats$err_reg)) # Voir si la variance est bien estimée
  
)

kable(Tab01)

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
#|n    |Methode    |Erreur de Monte Carlo                    |Autres                                   |
#|:----|:----------|:----------------------------------------|:----------------------------------------|
#|1000 |TNDDR_Mars |MCSE_bias          , 0.00496291544094027 |Bias               , -0.0305960863607402 |
#|-    |-          |MCSE_var           , 0.00191776654156114 |Var               , 0.0246305296739233   |
#|-    |-          |MCSE_mse           , 0.00190238093020336 |Mse               , 0.0255420196448433   |
#|-    |-          |%Cov , 0.897                             |Précision_var    , -4.66056457388497     |

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
