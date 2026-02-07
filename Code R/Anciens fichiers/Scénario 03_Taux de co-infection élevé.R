
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
# Scénario 03 (co-infection ~ 20% dans l'échantillon): Sous-scénario 01 : CV = 33%, 50% et 70% 
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
                        co_inf_para1 = 1, co_inf_para2 = -1, CV = j)
    
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
#|                 0.33|      0.4033667|    0.0107289|      0.4099743|    0.0107133|
#|                 0.50|      0.4061336|    0.0111115|      0.4127064|    0.0111236|
#|                 0.70|      0.4062554|    0.0084766|      0.4128896|    0.0084730|
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
    
    
    dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 1, co_inf_para2 = -1, CV = j, popsize = 1*10**6)
    
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
# nsim = 500
#|n    |Methode | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_err_reg| SD_coe_reg|
#|:----|:-------|--------------------:|---------:|---------:|------------:|----------:|
#|1000 |RegLog  |                 0.33| 0.4038068| 0.0676554|    0.1574603|  0.1668978|
#|-    |-       |                 0.50| 0.4070088| 0.0621061|    0.1439064|  0.1522959|
#|-    |-       |                 0.70| 0.3998882| 0.0665288|    0.1522310|  0.1655364|

# nsim = 1000
#|n    |Methode | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_err_reg| SD_coe_reg|
#|:----|:-------|--------------------:|---------:|---------:|------------:|----------:|
#|1000 |RegLog  |                 0.33| 0.4046778| 0.0693798|    0.1575425|  0.1710557|
#|-    |-       |                 0.50| 0.4053979| 0.0620926|    0.1438360|  0.1521965|
#|-    |-       |                 0.70| 0.4010034| 0.0673231|    0.1517962|  0.1680384|

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
# nsim = 500
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                                 |
#|:----|:--------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|
#|1000 |MCSE_bias           , 0.00302564214174351 , 0.00277747122473565 , 0.00297526011449047  |Bias               , 0.000440063580332162, 0.000875215197799761, -0.00636718266223973  |
#|-    |MCSE_var            , 0.000318984910006649, 0.000258543017493434, 0.000314403006555532 |Var                , 0.00457725518494713 , 0.00385717320211727 , 0.00442608637443892   |
#|-    |MCSE_mse            , 0.000319811692905452, 0.000259898477036569, 0.000304873832415989 |Mse                , 0.00456829433053197 , 0.00385022485735549 , 0.00445777521674437   |
#|-    |%Cov                , 0.942               , 0.938               , 0.93                 |Précision_var      , 0.00943750429553175 , 0.00838950265446697 , 0.0133053824344149    |

# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                                    |
#|:----|:--------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00219398298682414, 0.00196354070574119, 0.00212894302219092     |Bias                 , 0.00131107072176923  , -0.000735746811455107, -0.00525196343775264 |
#|-    |MCSE_var            , 0.000245274509253897, 0.000186710354250929, 0.00021880017339373  |Var                , 0.00481356134647376, 0.00385549210310263, 0.00453239839173543        |
#|-    |MCSE_mse            , 0.000247077425861482, 0.000185853036304505, 0.000213938749381439 |Mse                , 0.00481046669156477, 0.00385217793437009, 0.00455544911329518        |
#|-    |%Cov , 0.928, 0.931, 0.924                                                             |Précision_var      , 0.0135131686662445 , 0.00836051872330004, 0.0162421256685634         |
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
# nsim = 500
#|n    |Methode | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:-------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |IPW     |                 0.33| 0.4078760| 0.0703808|      0.1541582|  0.1705777|
#|-    |-       |                 0.50| 0.4111639| 0.0661214|      0.1471762|  0.1588096|
#|-    |-       |                 0.70| 0.4024455| 0.0682998|      0.1667137|  0.1692973|

# nsim = 1000
#|n    |Methode | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:-------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |IPW     |                 0.33| 0.4088095| 0.0719826|      0.1542225|  0.1749907|
#|-    |-       |                 0.50| 0.4093872| 0.0638516|      0.1472834|  0.1549365|
#|-    |-       |                 0.70| 0.4049040| 0.0717531|      0.1657223|  0.1780659|

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
# nsim = 500
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                                |
#|:----|:--------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------|
#|1000 |MCSE_bias           , 0.00314752726395533 , 0.00295703711391374 , 0.00305446039300825  |Bias               , -0.00209824140486708, -0.00154246137837266, -0.0104440701309048  |
#|-    |MCSE_var            , 0.000358858324179649, 0.000299199057437942, 0.000319330236413282 |Var                , 0.00495346393867108 , 0.00437203424653166 , 0.00466486414622807  |
#|-    |MCSE_mse            , 0.000354324397810899, 0.00029607551630344 , 0.000306456685336175 |Mse                , 0.00494795962778683 , 0.00436566936514237 , 0.00476461301883487  |
#|-    |%Cov                , 0.922               , 0.938               , 0.944                |Précision_var      , 0.0164195141005397  , 0.0116334123915284  , 0.00258359038717848  |

# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                                 |
#|:----|:--------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00227629004856144, 0.00201916624965861, 0.00226903228090763     |Bias                , -0.00116474361348917, -0.00331921026448057, -0.00798560251478048 |
#|-    |MCSE_var            , 0.000272328013114585, 0.000191297377943488, 0.000244622295509955 |Var                , 0.00518149638517986, 0.0040770323437604 , 0.0051485074918009      |
#|-    |MCSE_mse            , 0.000270621838921322, 0.000187607261557704, 0.000238006099334032 |Mse                , 0.00517767151647984, 0.00408397246819647, 0.00520712883183316     |
#|-    |%Cov , 0.911, 0.936, 0.94                                                              |Précision_var      , 0.0207681761596637 , 0.00765306903731389, 0.0123435459194287      |
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
# nsim = 5000
#|n    |Methode  | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_FR |                 0.33| 0.3986446| 0.0710786|      0.1738808|  0.1762666|
#|-    |-        |                 0.50| 0.4152428| 0.0709338|      0.1554130|  0.1670104|
#|-    |-        |                 0.70| 0.4265405| 0.0919418|      0.1882287|  0.2081234|

# nsim = 1000
#|n    |Methode  | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_FR |                 0.33| 0.3992712| 0.0735731|      0.1743335|  0.1826297|
#|-    |-        |                 0.50| 0.4126987| 0.0679682|      0.1536653|  0.1634474|
#|-    |-        |                 0.70| 0.4310103| 0.0994831|      0.2178156|  0.3436173|

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
# nsim = 500
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                             |
#|:----|:--------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias           , 0.00317873185660851 , 0.0031722566975061  , 0.00411176251346315  |Bias               , -0.0113296639756887, 0.00253638458100219, 0.0136508702305583  |
#|-    |MCSE_var            , 0.000393852044490638, 0.000469713997716929, 0.000731611922815384 |Var                , 0.0050521681081089 , 0.00503160627743614, 0.0084532954835604  |
#|-    |MCSE_mse            , 0.000373887507730239, 0.000476245225116873, 0.000793657276802743 |Mse                , 0.0051704250576947 , 0.00502797631162401, 0.00862273515064482 |
#|-    |%Cov                , 0.924               , 0.92                , 0.916                |Précision_var      , 0.00238576432866333, 0.0115974586419931 , 0.0198946570542484  |

# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                                  |
#|:----|:--------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00232658555787464, 0.00214934470563086, 0.00314593300719755     |Bias                 , -0.0107030902031696  , -7.67377945548597e-06, 0.0181207512058086 |
#|-    |MCSE_var            , 0.000301690920451303, 0.000222691334710049, 0.000833557266555452 |Var                , 0.00541300035811086, 0.00461968266362342, 0.00989689448577502      |
#|-    |MCSE_mse            , 0.000287890739419397, 0.000222680625842914, 0.000891276368136022 |Mse                , 0.00552214349764994, 0.00461506303984668, 0.0102153592155521       |
#|-    |%Cov , 0.917, 0.923, 0.929                                                             |Précision_var      , 0.0082961758629797 , 0.00978207885888094, 0.12580168161442         |
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
# nsim = 500
#|n    |Methode     | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:-----------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_Lasso |                 0.33| 0.4200205| 0.0736117|      0.0948602|  0.1741965|
#|-    |-           |                 0.50| 0.4199534| 0.0679559|      0.0829804|  0.1604950|
#|-    |-           |                 0.70| 0.4060410| 0.0859927|      0.0945273|  0.6234493|
  
# nsim = 1000
#|n    |Methode     | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:-----------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_Lasso |                 0.33| 0.4218521| 0.0770769|      0.0958852|  0.1795436|
#|-    |-           |                 0.50| 0.4185082| 0.0702191|      0.0840226|  0.2488750|
#|-    |-           |                 0.70| 0.4098520| 0.0899838|      0.0965747|  0.6262829|

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
# nsim = 500
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                 |Autres_(33%, 50%, 70%)                                                                 |
#|:----|:-------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|
#|1000 |MCSE_bias           , 0.00329201595070798 , 0.00303908075743517 , 0.00384571095194501 |Bias               , 0.0100462435377048  , 0.00724700890844254 , -0.00684861625662875  |
#|-    |MCSE_var            , 0.00037684581538166 , 0.000307722559543168, 0.00105162828356169 |Var                , 0.00541868450985789 , 0.00461800592510635 , 0.00739474636295492   |
#|-    |MCSE_mse            , 0.000403504816562827, 0.000324760370401283, 0.00105366516271297 |Mse                , 0.00550877415005705 , 0.00466128905137519 , 0.00742686041485957   |
#|-    |%Cov                , 0.726               , 0.69                , 0.65                |Précision_var      , 0.081406418500597   , 0.084030064721248   , 0.113596062035986     |
  
# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                             |
#|:----|:--------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00243738646221198, 0.00222052437899999, 0.00284553713056565     |Bias               , 0.0118778542163545 , 0.00580177224566053, -0.0030376018940147 |
#|-    |MCSE_var            , 0.000430552183220245, 0.000426163429639983, 0.000816738343135294 |Var                , 0.00594085276617422, 0.00493072851773328, 0.00809708156142778 |
#|-    |MCSE_mse            , 0.000455081895599498, 0.000433918445836256, 0.000816132655975758 |Mse                , 0.00607599533419302, 0.00495945835040606, 0.00809821150513287 |
#|-    |%Cov , 0.717, 0.69 , 0.636                                                             |Précision_var     , 0.0867444187315002, 0.0794247640290172, 0.24704261168542       |
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
# nsim = 500
#|n    |Methode    | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:----------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_Mars |                 0.33| 0.3946917| 0.1745081|      1.3068406|   1.443025|
#|-    |-          |                 0.50| 0.3866267| 0.1542198|      0.4273857|   1.515196|
#|-    |-          |                 0.70| 0.3989735| 0.1865545|      1.1004015|   1.654884|

# nsim = 1000
#|n    |Methode    | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:----------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_Mars |                 0.33| 0.3843854| 0.1707585|      0.8801430|   1.510912|
#|-    |-          |                 0.50| 0.3927772| 0.1490157|      0.4596358|   1.406908|
#|-    |-          |                 0.70| 0.4078505| 0.1792546|      0.5735527|   1.521206|

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
# nsim = 500
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                              |Autres_(33%, 50%, 70%)                                                            |
#----  |:----------------------------------------------------------------------------------|:---------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00780423983348915, 0.0068969208006454 , 0.00834297255850447 |Bias               , -0.015282584772547 , -0.0260796706150173, -0.013916087168636 |
#|-    |MCSE_var           , 0.00334838377446028, 0.00257244076948055, 0.00324111410625154 |Var                , 0.0304530796893094 , 0.0237837582651876 , 0.0348025955559793 |
#|-    |MCSE_mse           , 0.00327535087138251, 0.00261225221274743, 0.00321277162107997 |Mse                , 0.0306257309272609 , 0.024416339968045  , 0.0349266478469524 |
#|-    |%Cov               , 0.924              , 0.886              , 0.87                |Précision_var      , 0.136184078263234  , 1.08781037456674   , 0.554482875575702  |
  
# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                              |Autres_(33%, 50%, 70%)                                                                 |
#|:----|:----------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00539985879799289, 0.00471228897815777, 0.00566852931798278 |Bias                , -0.0255888612154266 , -0.019929203664652  , -0.00503910721745204 |
#|-    |MCSE_var           , 0.00221214502903715, 0.00181683901619541, 0.00214448848905477 |Var               , 0.0291584750382612, 0.0222056674136672, 0.0321322246288303         |
#|-    |MCSE_mse           , 0.00216145042793707, 0.00181803955238865, 0.00213788424740125 |Mse               , 0.0297841063815253, 0.0225806349049607, 0.0321254850057505         |
#|-    |%Cov , 0.914, 0.913, 0.866                                                         |Précision_var    , 0.63076935433168 , 0.947271893018713, 0.947653691134695             |
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
# nsim = 500
#|n    |Methode  | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_RN |                 0.33| 0.3820585| 0.1735947|      0.5927845|   1.584193|
#|-    |-        |                 0.50| 0.4072511| 0.1842058|      0.4845503|   1.674599|
#|-    |-        |                 0.70| 0.4127203| 0.2056823|      1.6716750|   1.870072|
  
# nsim = 1000
#|n    |Methode  | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_RN |                 0.33| 0.3812465| 0.1663434|      0.5246539|   1.511261|
#|-    |-        |                 0.50| 0.3942767| 0.1661963|      0.5070995|   1.662762|
#|-    |-        |                 0.70| 0.4211607| 0.2035735|      0.5567892|   1.857945|
  
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
# nsim = 500
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                              |Autres_(33%, 50%, 70%)                                                                     |
#|:----|:----------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00776339209385682, 0.00823793238284133, 0.00919839343847375 |Bias              , -0.0279157370995157  , -0.00545530622142681 , -0.000169286746305175    |
#|-    |MCSE_var           , 0.00289248848851221, 0.00310765817713159, 0.00342731720481951 |Var               , 0.0301351284014793   , 0.0339317649721329   , 0.0423052209244784       |
#|-    |MCSE_mse           , 0.002854046325862  , 0.00310015360637516, 0.0034269975450058  |Mse               , 0.0308541465224856   , 0.0338936618081582   , 0.0422206391406319       |
#|-    |%Cov               , 0.92               , 0.898              , 0.87                |Précision_var     , -0.416517907321245   , -0.317539874507732   , -1.46355160526974        |

# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                              |Autres_(33%, 50%, 70%)                                                             |
#|:----|:----------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00526024038667079, 0.00525558920793494, 0.00643755880682539 |Bias               , -0.0287278306253117, -0.0184296670859352, 0.00827111762182159 |
#|-    |MCSE_var           , 0.0019008663056102 , 0.00191034984858402, 0.00237199831069507 |Var               , 0.0276701289255625, 0.0276212179225622, 0.0414421633913351     |
#|-    |MCSE_mse           , 0.00188619337466732, 0.0019347585783995 , 0.00237756899387524 |Mse               , 0.0284677470490735, 0.027933249333538 , 0.0414691326146578     |
#|-    |%Cov , 0.904, 0.916, 0.869                                                         |Précision_var     , -0.342024290449187, -0.343652128516944, -0.213171870189708     |