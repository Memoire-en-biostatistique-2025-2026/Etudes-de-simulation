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
# Scénario 03 (co-infection ~ 40% dans l'échantillon): Sous-scénario 01 : CV = 33%, 50% et 70% 
#                                                                         I1_prev = 15% et 
#                                                                         I2_rev = 50%
################################################################################

# Listes des valeurs à tester (à discuter)

CV_liste <- c(0.33, 0.5, 0.7)

# Calcul des vraies valeurs 

set.seed(1) # Pour avoir toujours les mêmes germes

nsim1 <- 10 # Nombre de réplications

seeds_list <- sample(1:1000000, size = nsim1*3)

l_vraiRRc <- rep(NA, nsim1*3)
l_vraiRRm <- rep(NA, nsim1*3)

d <- 1
f <- nsim1

for(j in CV_liste) {
  
  for (i in d:f) {
    
    dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000,
                        co_inf_para1 = 1, co_inf_para2 = 2, CV = j)
    
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
  
  `Couverture vaccinale` = c(0.33, 0.5, 0.7),
  Mean_l_vraiRRc = c(mean(l_vraiRRc[1:nsim1]), mean(l_vraiRRc[(nsim1 + 1):(2*nsim1)]), mean(l_vraiRRc[(2*nsim1 + 1):(3*nsim1)])) ,
  SD_l_vraiRRc = c(sd(l_vraiRRc[1:nsim1]), sd(l_vraiRRc[(nsim1 + 1):(2*nsim1)]),sd(l_vraiRRc[(2*nsim1 + 1):(3*nsim1)])),
  Mean_l_vraiRRm =  c(mean(l_vraiRRm[1:nsim1]), mean(l_vraiRRm[(nsim1 + 1):(2*nsim1)]), mean(l_vraiRRm[(2*nsim1 + 1):(3*nsim1)])),
  SD_l_vraiRRm = c(sd(l_vraiRRm[1:nsim1]), sd(l_vraiRRm[(nsim1 + 1):(2*nsim1)]), sd(l_vraiRRm[(2*nsim1 + 1):(3*nsim1)]))
  
)

kable(Tab)
#| Couverture.vaccinale| Mean_l_vraiRRc| SD_l_vraiRRc| Mean_l_vraiRRm| SD_l_vraiRRm|
#|--------------------:|--------------:|------------:|--------------:|------------:|
#|                 0.33|      0.4123015|    0.0097757|      0.4188461|    0.0097158|
#|                 0.50|      0.4167544|    0.0061220|      0.4232299|    0.0061197|
#|                 0.70|      0.4157839|    0.0093774|      0.4222880|    0.0093950|

Tab_vrai4 <- Tab
l_vraiRRc_4 <- l_vraiRRc
l_vraiRRm_4 <- l_vraiRRm

################################################################################
# Initialiser des objets pour contenir les resultats

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

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

for(j in CV_liste) {
  
  for (i in d:f) {
    
    
    dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 1, co_inf_para2 = 2, CV = j, popsize = 1*10**6)
    
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

resultats$`Couverture vaccinale` <- c("0.33", rep("-", nsim - 1,), "0.5", rep("-", nsim - 1,), "0.7", rep("-", nsim - 1,))
resultats$`Prévalence-I1` <- c("0.15", rep("-", 3*nsim - 1,))
resultats$`Prévalence-I2` <- c("0.5", rep("-", 3*nsim - 1,))

resultats2$`Couverture vaccinale` <- c("0.33", rep("-", nsim - 1,), "0.5", rep("-", nsim - 1,), "0.7", rep("-", nsim - 1,))
resultats2$`Prévalence-I1` <- c("0.15", rep("-", 3*nsim - 1,))
resultats2$`Prévalence-I2` <- c("0.5", rep("-", 3*nsim - 1,))


resultats3$`Couverture vaccinale` <- c("0.33", rep("-", nsim - 1,), "0.5", rep("-", nsim - 1,), "0.7", rep("-", nsim - 1,))
resultats3$`Prévalence-I1` <- c("0.15", rep("-", 3*nsim - 1,))
resultats3$`Prévalence-I2` <- c("0.5", rep("-", 3*nsim - 1,))

################################################################################

########################## Régression logistique ###############################

##    - statistiques descriptives/valeur de couverture vaccinale

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("RegLog", "-", "-"),
  `Couverture vaccinale` = c(0.33, 0.5, 0.7),
  Mean_RRC = c(mean(resultats$RRc[1:nsim]), mean(resultats$RRc[(nsim + 1):(2*nsim)]), mean(resultats$RRc[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats$RRc[1:nsim]), sd(resultats$RRc[(nsim + 1):(2*nsim)]),sd(resultats$RRc[(2*nsim + 1):(3*nsim)])),
  Mean_err_reg =  c(mean(resultats$err_reg[1:nsim]), mean(resultats$err_reg[(nsim + 1):(2*nsim)]), mean(resultats$err_reg[(2*nsim + 1):(3*nsim)])),
  SD_coe_reg = c(sd(resultats$coe_reg[1:nsim]), sd(resultats$coe_reg[(nsim + 1):(2*nsim)]), sd(resultats$coe_reg[(2*nsim + 1):(3*nsim)]))
  
)

kable(Tab)
#|n    |Methode | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_err_reg| SD_coe_reg|
#|:----|:-------|--------------------:|---------:|---------:|------------:|----------:|
#|1000 |RegLog  |                 0.33| 0.4028611| 0.0680538|    0.1576657|  0.1675126|
#|-    |-       |                 0.50| 0.4049081| 0.0627918|    0.1444100|  0.1542629|
#|-    |-       |                 0.70| 0.4037218| 0.0645764|    0.1524579|  0.1594993|

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

Tab01$`Erreur de Monte Carlo_(33%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(33%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                                |
#|:----|:--------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00215204933194983, 0.00198565249208464, 0.00204208439465191     |Bias                , -0.00944043571395864, -0.0118463417265639 , -0.0120620949589071 |
#|-    |MCSE_var            , 0.000232134077724195, 0.000192485222175066, 0.000193357714973063 |Var                , 0.0046313163271457 , 0.00394281581932193, 0.00417010867488086    |
#|-    |MCSE_mse            , 0.000221911779933084, 0.000184190832584904, 0.000185214152713753 |Mse                , 0.00471580683728794, 0.00407920881580514, 0.00431143270100367    |
#|-    |%Cov , 0.925, 0.928, 0.935                                                             |Précision_var      , 0.00984685340128807, 0.00985287927112397, 0.00704136872143715    |
################################## IPW #########################################

##    - statistiques descriptives/valeur de couverture vaccinale

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("IPW", "-", "-"),
  `Couverture vaccinale` = c(0.33, 0.5, 0.7),
  Mean_RRC = c(mean(resultats2$RRm[1:nsim]), mean(resultats2$RRm[(nsim + 1):(2*nsim)]), mean(resultats2$RRm[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats2$RRm[1:nsim]), sd(resultats2$RRm[(nsim + 1):(2*nsim)]),sd(resultats2$RRm[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats2$var_log_RRm[1:nsim])), mean(sqrt(resultats2$var_log_RRm[(nsim + 1):(2*nsim)])), mean(sqrt(resultats2$var_log_RRm[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats2$RRm[1:nsim])), sd(log(resultats2$RRm[(nsim + 1):(2*nsim)])), sd(log(resultats2$RRm[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)
#|n    |Methode | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:-------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |IPW     |                 0.33| 0.4058422| 0.0713441|      0.1549935|  0.1729434|
#|-    |-       |                 0.50| 0.4087587| 0.0651012|      0.1488489|  0.1584082|
#|-    |-       |                 0.70| 0.4062747| 0.0669523|      0.1671205|  0.1633684|

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

Tab01$`Erreur de Monte Carlo_(33%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(33%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                                 |
#|:----|:--------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00225609929942411, 0.00205867932526031, 0.00211721826522642     |Bias               , -0.0130039188291619, -0.0144711966327258, -0.0160133192000249     |
#|-    |MCSE_var            , 0.000280161098577157, 0.000202674560759671, 0.000211092664069193 |Var                , 0.00508998404886195, 0.00423816056425426, 0.00448261318260837     |
#|-    |MCSE_mse            , 0.000263352700923345, 0.000193196990193264, 0.000198705543098631 |Mse                , 0.00525399596972852, 0.00444333793567302, 0.00473455696122765     |
#|-    |%Cov , 0.916, 0.934, 0.949                                                             |Précision_var       , 0.0179498488043712  , 0.00955936379137987 , -0.00375206385615698 |
################################## TNDDR #######################################

# Avec forêt aléatoire

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("TNDDR_FR", "-", "-"),
  `Couverture vaccinale` = c(0.33, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_RF[1:nsim]), mean(resultats3$RRm_RF[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_RF[1:nsim]), sd(resultats3$RRm_RF[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-RF`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-RF`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-RF`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_RF[1:nsim])), sd(log(resultats3$RRm_RF[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)
#|n    |Methode  | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_FR |                 0.33| 0.3963218| 0.0736876|      0.1737058|  0.1822741|
#|-    |-        |                 0.50| 0.4125983| 0.0685706|      0.1547204|  0.1639618|
#|-    |-        |                 0.70| 0.4311832| 0.0915385|      0.1946651|  0.2756366|

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

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

Tab01$`Erreur de Monte Carlo_(33%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(33%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                             |
#|:----|:--------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00233020522906028, 0.00216839129440181, 0.00289470213176788     |Bias               , -0.0225243150724743, -0.0106315274831344, 0.00889517488204461 |
#|-    |MCSE_var            , 0.000293102293309903, 0.000252955636962486, 0.000726369293706061 |Var                , 0.00542985640953986, 0.00470192080563756, 0.00837930043166148 |
#|-    |MCSE_mse            , 0.000267190132524766, 0.00024031101840862 , 0.00075113692550728  |Mse                , 0.00593177132261441, 0.00481024826145656, 0.00845004526741198 |
#|-    |%Cov , 0.886, 0.9  , 0.931                                                             |Précision_var      , 0.00856827195819557, 0.00924134570813803, 0.0809714977946986  |
################################################################################
# Régression Lasso

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("TNDDR_Lasso", "-", "-"),
  `Couverture vaccinale` = c(0.33, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_Lasso[1:nsim]), mean(resultats3$RRm_Lasso[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_Lasso[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_Lasso[1:nsim]), sd(resultats3$RRm_Lasso[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_Lasso[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-Lasso`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-Lasso`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-Lasso`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_Lasso[1:nsim])), sd(log(resultats3$RRm_Lasso[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_Lasso[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)
#|n    |Methode     | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:-----------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_Lasso |                 0.33| 0.4196953| 0.0793613|      0.0971627|  0.2616696|
#|-    |-           |                 0.50| 0.4190336| 0.0747945|      0.0872568|  0.1679833|
#|-    |-           |                 0.70| 0.4123291| 0.0747853|      0.0857517|  0.3187495|

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

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

Tab01$`Erreur de Monte Carlo_(33%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(33%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                                 |
#|:----|:--------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00250962605388004, 0.00236520904171117, 0.00236491980388633     |Bias                , 0.000849182365890155, -0.00419623668760694, -0.00995897027595238 |
#|-    |MCSE_var            , 0.000561973689786079, 0.000615813946684188, 0.000472627479571183 |Var                , 0.0062982229303135 , 0.00559421381099227, 0.00559284567881374     |
#|-    |MCSE_mse            , 0.000563528693176848, 0.000606273302753374, 0.000465398454143283 |Mse                , 0.00629264581807372, 0.00560622799951969, 0.00568643392209223     |
#|-    |%Cov , 0.7  , 0.678, 0.646                                                             |Précision_var     , 0.0851113569427504, 0.0767049242407531, 0.189884899020962          |
################################################################################
# Régression avec splines: MARS

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("TNDDR_Mars", "-", "-"),
  `Couverture vaccinale` = c(0.33, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_Mars[1:nsim]), mean(resultats3$RRm_Mars[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_Mars[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_Mars[1:nsim]), sd(resultats3$RRm_Mars[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_Mars[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-Mars`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-Mars`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-Mars`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_Mars[1:nsim])), sd(log(resultats3$RRm_Mars[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_Mars[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)
#|n    |Methode    | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:----------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_Mars |                 0.33| 0.3870082| 0.1730713|      0.5865824|   1.367202|
#|-    |-          |                 0.50| 0.3935860| 0.1536511|      0.4177895|   1.370230|
#|-    |-          |                 0.70| 0.4087434| 0.1741033|      0.6435008|   1.537561|

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

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

Tab01$`Erreur de Monte Carlo_(33%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(33%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                              |Autres_(33%, 50%, 70%)                                                             |
#|:----|:----------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.0054729945736897 , 0.00485887506355917, 0.00550563121134277 |Bias               , -0.0318379517364126, -0.0296438975319593, -0.0135446200625786 |
#|-    |MCSE_var           , 0.0023418641494646 , 0.00184834641231137, 0.0021470282621766  |Var               , 0.0299536696036369, 0.0236086668832772, 0.0303119750353117     |
#|-    |MCSE_mse           , 0.00223252437439856, 0.00185127085628192, 0.0021331514108336  |Mse               , 0.0309373711048034, 0.0244638188772792, 0.030465119792916      |
#|-    |%Cov , 0.885, 0.888, 0.879                                                         |Précision_var    , 0.780619730474507, 0.952440911653342, 0.894060407698121         |
################################################################################
# Réseaux de neurones

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("TNDDR_RN", "-", "-"),
  `Couverture vaccinale` = c(0.33, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_RN[1:nsim]), mean(resultats3$RRm_RN[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_RN[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_RN[1:nsim]), sd(resultats3$RRm_RN[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_RN[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-RN`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-RN`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-RN`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_RN[1:nsim])), sd(log(resultats3$RRm_RN[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_RN[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)
#|n    |Methode  | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_RN |                 0.33| 0.3679109| 0.1767599|      0.6463732|   1.755971|
#|-    |-        |                 0.50| 0.3962882| 0.1667017|      0.5219108|   1.623689|
#|-    |-        |                 0.70| 0.4260689| 0.2104465|      0.5759940|   1.884714|

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

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

Tab01$`Erreur de Monte Carlo_(33%, 50%, 70%)` = list(
  
  list(name = "MCSE_bias", value = MCSE_biais),
  list(name = "MCSE_var", value = MCSE_var),
  list(name = "MCSE_mse", value = MCSE_mse),
  list(name = "%Cov", value = coverage)
  
)

Tab01$`Autres_(33%, 50%, 70%)` = list(
  
  list(name = "Bias", value = Biais),
  list(name = "Var", value = Var),
  list(name = "Mse", value = Mse),
  list(name = "Précision_var", value = Precision_var) # Voir si la variance est bien estimée
  
)

kable(Tab01)
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                              |Autres_(33%, 50%, 70%)                                                             |
#|:----|:----------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00558963758886933, 0.00527157089400848, 0.00665490262191544 |Bias               , -0.0509352130049855, -0.0269416719589455, 0.00378080958478005 |
#|-    |MCSE_var           , 0.00205852129722422, 0.00188272013037291, 0.00248953005087445 |Var               , 0.0312440483749009, 0.0277894596905573, 0.0442877289071769     |
#|-    |MCSE_mse           , 0.00207460874973899, 0.00191814373285877, 0.00249329351212924 |Mse               , 0.0338072002503893, 0.0284875239188102, 0.0442577356993861     |
#|-    |%Cov , 0.894, 0.882, 0.875                                                         |Précision_var     , -0.464099138323349, -0.357949029241477, -0.30035744644519      |