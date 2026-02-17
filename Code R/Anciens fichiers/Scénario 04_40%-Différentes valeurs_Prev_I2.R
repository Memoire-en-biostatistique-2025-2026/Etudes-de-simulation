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

#| Prévalence.de.I2| Mean_l_vraiRRc| SD_l_vraiRRc| Mean_l_vraiRRm| SD_l_vraiRRm|
#|----------------:|--------------:|------------:|--------------:|------------:|
#|              0.3|      0.4076055|    0.0081666|      0.4117046|    0.0081696|
#|              0.5|      0.4145370|    0.0077920|      0.4210375|    0.0078156|
#|              0.7|      0.4121933|    0.0066092|      0.4212769|    0.0065928|

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

for(j in I2_prev_liste) {
  
  for (i in d:f) {
    
    
    dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 1, co_inf_para2 = 2, I2_prev = j, popsize = 1*10**6)
    
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
  
  n = c("1000", "-", "-"), Methode = c("RegLog", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats$RRc[1:nsim]), mean(resultats$RRc[(nsim + 1):(2*nsim)]), mean(resultats$RRc[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats$RRc[1:nsim]), sd(resultats$RRc[(nsim + 1):(2*nsim)]),sd(resultats$RRc[(2*nsim + 1):(3*nsim)])),
  Mean_err_reg =  c(mean(resultats$err_reg[1:nsim]), mean(resultats$err_reg[(nsim + 1):(2*nsim)]), mean(resultats$err_reg[(2*nsim + 1):(3*nsim)])),
  SD_coe_reg = c(sd(resultats$coe_reg[1:nsim]), sd(resultats$coe_reg[(nsim + 1):(2*nsim)]), sd(resultats$coe_reg[(2*nsim + 1):(3*nsim)]))
  
)

kable(Tab)

#|n    |Methode | Prévalence.I2|  Mean_RRC|    SD_RRc| Mean_err_reg| SD_coe_reg|
#|:----|:-------|-------------:|---------:|---------:|------------:|----------:|
#|1000 |RegLog  |           0.3| 0.4048538| 0.0658187|    0.1513492|  0.1621767|
#|-    |-       |           0.5| 0.4062589| 0.0669741|    0.1575880|  0.1659333|
#|-    |-       |           0.7| 0.4063858| 0.0695150|    0.1673344|  0.1711335|

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
#|n    |Erreur de Monte Carlo_(30%, 50%, 70%)                                                  |Autres_(30%, 50%, 70%)                                                                 |
#|:----|:--------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00208136918406825, 0.0021179076263142 , 0.00219825842389896     |Bias                , -0.00275161742070568, -0.00827813995447124, -0.00580750371123345 |
#|-    |MCSE_var            , 0.000208987547762832, 0.000207473217375701, 0.00021413751085504  |Var                , 0.00433209768038895, 0.00448553271359983, 0.00483234009824276     |
#|-    |MCSE_mse            , 0.000205971472395468, 0.000202403261489562, 0.000208568698625228 |Mse                , 0.00433533698113849, 0.00454957478199204, 0.00486123485750051     |
#|-    |%Cov , 0.931, 0.933, 0.951                                                             |Précision_var      , 0.010827414153898  , 0.00834527864378354, 0.00379911662696955     |
################################## IPW #########################################

##    - statistiques descriptives/valeur de couverture vaccinale

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("IPW", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats2$RRm[1:nsim]), mean(resultats2$RRm[(nsim + 1):(2*nsim)]), mean(resultats2$RRm[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats2$RRm[1:nsim]), sd(resultats2$RRm[(nsim + 1):(2*nsim)]),sd(resultats2$RRm[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats2$var_log_RRm[1:nsim])), mean(sqrt(resultats2$var_log_RRm[(nsim + 1):(2*nsim)])), mean(sqrt(resultats2$var_log_RRm[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats2$RRm[1:nsim])), sd(log(resultats2$RRm[(nsim + 1):(2*nsim)])), sd(log(resultats2$RRm[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)

#|n    |Methode | Prévalence.I2|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:-------|-------------:|---------:|---------:|--------------:|----------:|
#|1000 |IPW     |           0.3| 0.4076247| 0.0661941|      0.1405213|  0.1620890|
#|-    |-       |           0.5| 0.4093138| 0.0682268|      0.1547712|  0.1676519|
#|-    |-       |           0.7| 0.4086417| 0.0710454|      0.1723634|  0.1744274|

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

#|n    |Erreur de Monte Carlo_(30%, 50%, 70%)                                                  |Autres_(30%, 50%, 70%)                                                                |
#|:----|:--------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.0020932423060096 , 0.0021575196932256 , 0.00224665321589655     |Bias                , -0.00407984739630424, -0.0117237095483646 , -0.0126351832083589 |
#|-    |MCSE_var            , 0.000213560024032281, 0.000212883269875989, 0.0002173693830592   |Var                , 0.0043816633516684 , 0.00465489122665631, 0.00504745067249832    |
#|-    |MCSE_mse            , 0.000209367604898615, 0.000206974202781009, 0.000209846834082756 |Mse                , 0.00439392684309387, 0.00478768170100407, 0.00520205107653461    |
#|-    |%Cov , 0.915, 0.926, 0.952                                                             |Précision_var      , 0.0215676963110846 , 0.0128806365261095 , 0.00206391462632227    |
################################## TNDDR #######################################

# Avec forêt aléatoire

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("TNDDR_FR", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_RF[1:nsim]), mean(resultats3$RRm_RF[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_RF[1:nsim]), sd(resultats3$RRm_RF[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-RF`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-RF`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-RF`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_RF[1:nsim])), sd(log(resultats3$RRm_RF[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_RF[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)

#|n    |Methode  | Prévalence.I2|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|-------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_FR |           0.3| 0.3992593| 0.0699942|      0.1635944|  0.1752002|
#|-    |-        |           0.5| 0.3976864| 0.0708148|      0.1748429|  0.1794808|
#|-    |-        |           0.7| 0.3983212| 0.0752239|      0.1885119|  0.1893352|

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
#|n    |Erreur de Monte Carlo_(30%, 50%, 70%)                                                  |Autres_(30%, 50%, 70%)                                                                 |
#|:----|:--------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00221341033824346, 0.00223936074463963, 0.00237878958321613     |Bias               , -0.0124452418572777, -0.0233511285178349, -0.0229557359465375     |
#|-    |MCSE_var            , 0.000244420927094726, 0.000232993360037398, 0.000266589146040026 |Var                , 0.00489918532544304, 0.00501473654463296, 0.00565863988121756     |
#|-    |MCSE_mse            , 0.000234170485265446, 0.00023088002425401 , 0.000256733690962355 |Mse                , 0.00504917018500373, 0.00555499701114477, 0.0061799470541835      |
#|-    |%Cov , 0.905, 0.891, 0.903                                                             |Précision_var       , 0.011605765137861   , 0.00463788409647725 , 0.000823298419201896 |
################################################################################
# Régression Lasso

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("TNDDR_Lasso", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_Lasso[1:nsim]), mean(resultats3$RRm_Lasso[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_Lasso[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_Lasso[1:nsim]), sd(resultats3$RRm_Lasso[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_Lasso[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-Lasso`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-Lasso`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-Lasso`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_Lasso[1:nsim])), sd(log(resultats3$RRm_Lasso[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_Lasso[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)
#|n    |Methode     | Prévalence.I2|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:-----------|-------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_Lasso |           0.3| 0.4184150| 0.0692994|      0.1035893|  0.1646217|
#|-    |-           |           0.5| 0.4229011| 0.0716914|      0.0945311|  0.1705436|
#|-    |-           |           0.7| 0.4323054| 0.1159272|      1.7408578|  0.7163666|

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
#|n    |Erreur de Monte Carlo_(30%, 50%, 70%)                                                 |Autres_(30%, 50%, 70%)                                                            |
#|:----|:-------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00219144098898361, 0.0022670807530425 , 0.0036659412647822     |Bias               , 0.00671039907988275, 0.00186360011872616, 0.0110284298355689 |
#|-    |MCSE_var            , 0.000239494834845326, 0.000240908063904668, 0.00145428483039203 |Var                , 0.00480241360819747, 0.00513965514081576, 0.0134391253568329 |
#|-    |MCSE_mse            , 0.000250311979802593, 0.00024307353888934 , 0.00148853022423214 |Mse                , 0.00484264065040057, 0.00513798849107746, 0.0135473124961141 |
#|-    |%Cov , 0.78 , 0.73 , 0.686                                                            |Précision_var     , 0.0716109148881023, 0.0849497049760303, -1.55152251116834     |
################################################################################
# Régression avec splines: MARS

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("TNDDR_Mars", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_Mars[1:nsim]), mean(resultats3$RRm_Mars[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_Mars[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_Mars[1:nsim]), sd(resultats3$RRm_Mars[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_Mars[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-Mars`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-Mars`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-Mars`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_Mars[1:nsim])), sd(log(resultats3$RRm_Mars[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_Mars[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)
#|n    |Methode    | Prévalence.I2|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:----------|-------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_Mars |           0.3| 0.3922834| 0.1587123|      0.8066186|   1.231961|
#|-    |-          |           0.5| 0.3836907| 0.1617596|      0.6751800|   1.436712|
#|-    |-          |           0.7| 0.3723224| 0.1727830|      0.5517511|   1.587666|
  

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
#|n    |Erreur de Monte Carlo_(30%, 50%, 70%)                                              |Autres_(30%, 50%, 70%)                                                             |
#|:----|:----------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00501892248111425, 0.00511528697355722, 0.00546387750570997 |Bias               , -0.0194211432889249, -0.0373468891796407, -0.0489545194389134 |
#|-    |MCSE_var           , 0.00205748951351715, 0.00200164232210413, 0.00213871907740257 |Var               , 0.025189582871434 , 0.0261661608218441, 0.0298539573974034     |
#|-    |MCSE_mse           , 0.00200579399172471, 0.00196762484110307, 0.002093327965715   |Mse               , 0.0255415740952115, 0.0275347847924187, 0.032220648413501      |
#|-    |%Cov , 0.912, 0.896, 0.899                                                         |Précision_var    , 0.425342263110731, 0.761531922803584, 1.03591528159927          |
################################################################################
# Réseaux de neurones

Tab <- data.frame( 
  
  n = c("1000", "-", "-"), Methode = c("TNDDR_RN", "-", "-"),
  `Prévalence-I2` = c(0.3, 0.5, 0.7),
  Mean_RRC = c(mean(resultats3$RRm_RN[1:nsim]), mean(resultats3$RRm_RN[(nsim + 1):(2*nsim)]), mean(resultats3$RRm_RN[(2*nsim + 1):(3*nsim)])) ,
  SD_RRc = c(sd(resultats3$RRm_RN[1:nsim]), sd(resultats3$RRm_RN[(nsim + 1):(2*nsim)]),sd(resultats3$RRm_RN[(2*nsim + 1):(3*nsim)])),
  Mean_sd_LOgRRm =  c(mean(sqrt(resultats3$`var_log_RRm-RN`[1:nsim])), mean(sqrt(resultats3$`var_log_RRm-RN`[(nsim + 1):(2*nsim)])), mean(sqrt(resultats3$`var_log_RRm-RN`[(2*nsim + 1):(3*nsim)]))),
  SD_log_RRm = c(sd(log(resultats3$RRm_RN[1:nsim])), sd(log(resultats3$RRm_RN[(nsim + 1):(2*nsim)])), sd(log(resultats3$RRm_RN[(2*nsim + 1):(3*nsim)])))
  
)

kable(Tab)
#|n    |Methode  | Prévalence.I2|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|-------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_RN |           0.3| 0.3869586| 0.1300839|      0.3816416|   1.143290|
#|-    |-        |           0.5| 0.3717350| 0.1637631|      0.6984780|   1.616861|
#|-    |-        |           0.7| 0.3592959| 0.1939743|      1.5777823|   1.945314|

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
#|n    |Erreur de Monte Carlo_(30%, 50%, 70%)                                              |Autres_(30%, 50%, 70%)                                                             |
#|:----|:----------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00411361379344619, 0.00517864476825501, 0.00613400716900349 |Bias               , -0.0247459332273098, -0.0493025173545057, -0.0619810553617557 |
#|-    |MCSE_var           , 0.0015140330258995 , 0.00174661894187441, 0.00218620413921531 |Var               , 0.0169218184416307, 0.026818361635775 , 0.0376260439493863     |
#|-    |MCSE_mse           , 0.00151208757917455, 0.00182968307849597, 0.00222010740333567 |Mse               , 0.0175172578344796, 0.0292222814916305, 0.0414300691291939     |
#|-    |%Cov , 0.919, 0.883, 0.889                                                         |Précision_var     , -0.206441441919406, -0.5189971929018  , -1.3884471029762       |
