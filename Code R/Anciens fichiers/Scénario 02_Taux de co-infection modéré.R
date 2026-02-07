
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
# Scénario 02 (co-infection ~ 10% dans l'échantillon): Sous-scénario 01 : CV = 33%, 50% et 70% 
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
                      co_inf_para1 = 2, co_inf_para2 = -2, CV = j)
  
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

Tab_Vrai <- data.frame( 
  
  `Couverture vaccinale` = c(0.33, 0.5, 0.7),
  Mean_l_vraiRRc = c(mean(l_vraiRRc[1:nsim1]), mean(l_vraiRRc[(nsim1 + 1):(2*nsim1)]), mean(l_vraiRRc[(2*nsim1 + 1):(3*nsim1)])) ,
  SD_l_vraiRRc = c(sd(l_vraiRRc[1:nsim1]), sd(l_vraiRRc[(nsim1 + 1):(2*nsim1)]),sd(l_vraiRRc[(2*nsim1 + 1):(3*nsim1)])),
  Mean_l_vraiRRm =  c(mean(l_vraiRRm[1:nsim1]), mean(l_vraiRRm[(nsim1 + 1):(2*nsim1)]), mean(l_vraiRRm[(2*nsim1 + 1):(3*nsim1)])),
  SD_l_vraiRRm = c(sd(l_vraiRRm[1:nsim1]), sd(l_vraiRRm[(nsim1 + 1):(2*nsim1)]), sd(l_vraiRRm[(2*nsim1 + 1):(3*nsim1)]))
  
)

kable(Tab_Vrai)

#| Couverture.vaccinale| Mean_l_vraiRRc| SD_l_vraiRRc| Mean_l_vraiRRm| SD_l_vraiRRm|
#|--------------------:|--------------:|------------:|--------------:|------------:|
#|                 0.33|      0.4116813|    0.0076165|      0.4180096|    0.0076286|
#|                 0.50|      0.4131110|    0.0111406|      0.4194260|    0.0111263|
#|                 0.70|      0.4143477|    0.0104384|      0.4206930|    0.0104321|

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
    
  
    dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 2, co_inf_para2 = -2, CV = j, popsize = 1*10**6)
  
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
#|1000 |RegLog  |                 0.33| 0.4250327| 0.0696351|    0.1565009|  0.1632441|
#|-    |-       |                 0.50| 0.4169054| 0.0615427|    0.1429681|  0.1466785|
#|-    |-       |                 0.70| 0.4210229| 0.0681662|    0.1508180|  0.1627938|

# nsim = 1000
#|n    |Methode | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_err_reg| SD_coe_reg|
#|:----|:-------|--------------------:|---------:|---------:|------------:|----------:|
#|1000 |RegLog  |                 0.33| 0.4230303| 0.0685047|    0.1563896|  0.1620452|
#|-    |-       |                 0.50| 0.4141558| 0.0626125|    0.1430297|  0.1511648|
#|-    |-       |                 0.70| 0.4161979| 0.0688661|    0.1508837|  0.1656081|

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
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                             |
#|:----|:--------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias           , 0.00311417455701115 , 0.0027522727612552  , 0.00304848435162244  |Bias               , 0.013351380490392  , 0.00379441865524693, 0.00667521202457666 |
#|-    |MCSE_var            , 0.000359514377450476, 0.000275031648173041, 0.00028651440389443  |Var                , 0.00484904158576781, 0.00378750267617366, 0.00464662842104345 |
#|-    |MCSE_mse            , 0.000392656213308733, 0.000282029975280511, 0.000298258029974906 |Mse                , 0.0050176028635955 , 0.0037943252837526 , 0.00468189361977442 |
#|-    |%Cov                , 0.944               , 0.928               , 0.934                |Précision_var      , 0.00674322617009104, 0.00371037544669611, 0.0119757823668358  |

# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                             |
#|:----|:--------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00216630980391496, 0.00197998174497109, 0.00217773617106415     |Bias               , 0.0113489504943576 , 0.00104480538496637, 0.00185021893957366 |
#|-    |MCSE_var            , 0.000226306353449111, 0.000174425151055008, 0.000237041126637361 |Var                , 0.00469289816653806, 0.00392032771041875, 0.00474253483076113 |
#|-    |MCSE_mse            , 0.000243947963307224, 0.000175523490949166, 0.000239391801070651 |Mse                , 0.0048170039456949 , 0.00391749900100078, 0.00474121560605472 |
#|-    |%Cov , 0.937, 0.939, 0.925                                                             |Précision_var      , 0.00565559999623541, 0.00813506457029864, 0.0147243481948149  |
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
#|1000 |IPW     |                 0.33| 0.4300749| 0.0726145|      0.1519213|  0.1680536|
#|-    |-       |                 0.50| 0.4218344| 0.0639344|      0.1445933|  0.1502887|
#|-    |-       |                 0.70| 0.4252008| 0.0715944|      0.1617404|  0.1706173|

# nsim = 1000

help("calc_absolute") # calculer les différentes mesures de performance
#|n    |Methode | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:-------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |IPW     |                 0.33| 0.4267937| 0.0707502|      0.1518352|  0.1653077|
#|-    |-       |                 0.50| 0.4171922| 0.0643166|      0.1450969|  0.1544051|
#|-    |-       |                 0.70| 0.4204609| 0.0728564|      0.1623887|  0.1740626|
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
# nsim =500
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                             |
#|:----|:--------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias           , 0.00324741913753016 , 0.00285923299128673 , 0.00320179753718091  |Bias               , 0.0120652988263494 , 0.00240838918096647, 0.00450778088755277 |
#|-    |MCSE_var            , 0.000417375084519998, 0.000306690569254914, 0.000317385307111493 |Var                , 0.00527286552739855, 0.00408760664923124, 0.00512575373454888 |
#|-    |MCSE_mse            , 0.000448079748970355, 0.000311334600598532, 0.000323753244639615 |Mse                , 0.00540789123211286, 0.00408523177437977, 0.00513582231560997 |
#|-    |%Cov                , 0.93                , 0.948               , 0.938                |Précision_var      , 0.0161323057523936 , 0.00569548461175934, 0.00887681522278336 |

# nsim =1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                                     |
#|:----|:--------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00223731892363462, 0.00203386869798314, 0.00230392172957643     |Bias                 , 0.00878408879013964  , -0.00223376798834241 , -0.000232050306398779 |
#|-    |MCSE_var            , 0.00025740035534772 , 0.000184213800820116, 0.000262689273502226 |Var                , 0.00500559596605356, 0.00413662188063561, 0.00530805533601443         |
#|-    |MCSE_mse            , 0.000271501658478269, 0.000182217633250731, 0.000262404630910948 |Mse                , 0.00507775058596056, 0.00413747497818072, 0.00530280112802312         |
#|-    |%Cov , 0.932, 0.937, 0.934                                                             |Précision_var      , 0.0134725070817817 , 0.00930816819712577, 0.0116739358728295        
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
# nsim = 500
#|n    |Methode  | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_FR |                 0.33| 0.4198753| 0.0773074|      0.1740871|  0.1823423|
#|-    |-        |                 0.50| 0.4245751| 0.0648730|      0.1520627|  0.1524510|
#|-    |-        |                 0.70| 0.4500305| 0.0954087|      0.1870237|  0.2018623|
  
# nsim = 1000
#|n    |Methode  | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_FR |                 0.33| 0.4155387| 0.0744678|      0.1734711|  0.1780439|
#|-    |-        |                 0.50| 0.4215293| 0.0687211|      0.1525260|  0.1633625|
#|-    |-        |                 0.70| 0.4458225| 0.0951370|      0.1870380|  0.2037870|

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
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                 |Autres_(33%, 50%, 70%)                                                               |
#|:----|:-------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------|
#|1000 |MCSE_bias           , 0.00345729323153052 , 0.00290120962926638 , 0.00426680881737249 |Bias               , 0.00186573984695027, 0.00514913155843261 , 0.0293375374869498   |
#|-    |MCSE_var            , 0.00051242062296961 , 0.000294922815895986, 0.00108977663136514 |Var                , 0.00597643824439338, 0.00420850865647398 , 0.00910282874200379  |
#|-    |MCSE_mse            , 0.000517147268648621, 0.000304522652524815, 0.00123583529151197 |Mse                , 0.00596796635308109, 0.00422660519496708 , 0.00994531419031797  |
#|-    |%Cov                , 0.936               , 0.928               , 0.95                |Précision_var      , 0.00825521485448824, 0.000388280672584451, 0.0148385806456601   |
  
# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                               |
#|:----|:--------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00235487893807734, 0.00217315202360962, 0.00300849590457903     |Bias                , -0.00247091305638913, 0.00210332249762646 , 0.0251295504913137 |
#|-    |MCSE_var            , 0.000310312815425886, 0.000217755781846677, 0.000689492173208832 |Var                , 0.00554545481300025, 0.00472258971771861, 0.00905104760786877   |
#|-    |MCSE_mse            , 0.000306660681852217, 0.000220313321273955, 0.00077803088706692  |Mse                , 0.00554601476951948, 0.00472229109352991, 0.00967349086815639   |
#|-    |%Cov , 0.929, 0.926, 0.934                                                             |Précision_var      , 0.00457284460042992, 0.0108364736102791 , 0.0167489454205554    |
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
#|1000 |TNDDR_Lasso |                 0.33| 0.4434795| 0.0753724|      0.0947451|  0.1694810|
#|-    |-           |                 0.50| 0.4323497| 0.0705614|      0.0851761|  0.1570158|
#|-    |-           |                 0.70| 0.4298705| 0.0777833|      0.0871757|  0.4194302|
  
# nsim = 1000
#|n    |Methode     | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:-----------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_Lasso |                 0.33| 0.4398381| 0.0736115|      0.0948302|  0.1672303|
#|-    |-           |                 0.50| 0.4269491| 0.0714323|      0.0848438|  0.1616012|
#|-    |-           |                 0.70| 0.4264564| 0.0854813|      0.0914374|  0.4629420|
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
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                             |
#|:----|:--------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias           , 0.00337075643716001 , 0.00315560248004053 , 0.00347857615780419  |Bias               , 0.0254698996814597 , 0.0129236922962461 , 0.00917754525032521 |
#|-    |MCSE_var            , 0.000412658863956251, 0.000702958774133771, 0.000600978346807581 |Var                , 0.00568099947932782, 0.00497891350601897, 0.00605024604282189 |
#|-    |MCSE_mse            , 0.000492699334485432, 0.00074147893721484 , 0.000590798600292211 |Mse                , 0.00631835327015279, 0.00513597750157498, 0.00612237288755801 |
#|-    |%Cov                , 0.746               , 0.728               , 0.664                |Précision_var      , 0.0875971980647791 , 0.0672749622587076 , 0.114686626050516   |
  
# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                                  |Autres_(33%, 50%, 70%)                                                             |
#|:----|:--------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00232780080707228, 0.00225888661754905, 0.00270315525527162     |Bias               , 0.0218284780716697 , 0.00752314470819754, 0.00576342221368131 |
#|-    |MCSE_var            , 0.000266593902468301, 0.000497721047928933, 0.000731441723466429 |Var                , 0.00541865659740636, 0.0051025687509422 , 0.00730704833410257 |
#|-    |MCSE_mse            , 0.000311054894787999, 0.000512649744793659, 0.000737476397362216 |Mse                , 0.00588972039573432, 0.00515406388849174, 0.00733295832138163 |
#|-    |%Cov , 0.744, 0.705, 0.665                                                             |Précision_var     , 0.0832136937374597, 0.0785187585132607, 0.11234959814632       |
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
#|1000 |TNDDR_Mars |                 0.33| 0.3941600| 0.1669337|      0.6128052|   1.535539|
#|-    |-          |                 0.50| 0.4103104| 0.1550793|      0.3909655|   1.332055|
#|-    |-          |                 0.70| 0.4293446| 0.1719070|      0.4679842|   1.401758|
  
# nsim = 1000
#|n    |Methode    | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:----------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_Mars |                 0.33| 0.3980158| 0.1752132|      0.5918660|   1.558413|
#|-    |-          |                 0.50| 0.4049988| 0.1550671|      0.5455736|   1.451525|
#|-    |-          |                 0.70| 0.4199309| 0.1721314|      0.5655668|   1.497915|

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
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                              |Autres_(33%, 50%, 70%)                                                                |
#|:----|:----------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00746550214269901, 0.00693535492990164, 0.00768791307027893 |Bias              , -0.0238496142561721 , -0.00911560066542072, 0.00865160874681725   |
#|-    |MCSE_var           , 0.00285731409307373, 0.00277850796670404, 0.00290655062233337 |Var               , 0.0278668611213218  , 0.0240495740018555  , 0.0295520036880828    |
#|-    |MCSE_mse           , 0.00286367622884166, 0.00276572854575415, 0.00291896796453859 |Mse               , 0.0283799314992473  , 0.0240845690293432  , 0.0295677500146146    |
#|-    |%Cov               , 0.926              , 0.894              , 0.918               |Précision_var     , 0.922734095980471   , 0.941089057731959   , 0.933774327268879     |
  
# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                              |Autres_(33%, 50%, 70%)                                                                     |
#|:----|:----------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00554072830225283, 0.00490365100901632, 0.00544327314250841 |Bias                 , -0.0199937695491645  , -0.0144272153882096  , -0.000762098840571546 |
#|-    |MCSE_var           , 0.00222965780971019, 0.00195586460485949, 0.00201904075672403 |Var               , 0.0306996701193855, 0.0240457932182268, 0.0296292225039533             |
#|-    |MCSE_mse           , 0.00220668500441829, 0.00195110067017057, 0.00201896406380439 |Mse               , 0.0310687212700512, 0.0242298919688664, 0.0296001740760922             |
#|-    |%Cov , 0.929, 0.915, 0.884                                                         |Précision_var    , 0.966547559900786, 0.905950919188625, 0.932347946697915                 |
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
#|1000 |TNDDR_RN |                 0.33| 0.3892476| 0.1719195|    176.2534003|   1.721679|
#|-    |-        |                 0.50| 0.4169466| 0.1668164|      0.3355555|   1.510608|
#|-    |-        |                 0.70| 0.4401697| 0.2121747|      0.7360616|   1.854088|
  
# nsim = 1000
#|n    |Methode  | Couverture.vaccinale|  Mean_RRC|    SD_RRc| Mean_sd_LOgRRm| SD_log_RRm|
#|:----|:--------|--------------------:|---------:|---------:|--------------:|----------:|
#|1000 |TNDDR_RN |                 0.33| 0.3871027| 0.1686434|     88.5385842|   1.654101|
#|-    |-        |                 0.50| 0.4062819| 0.1630770|      0.5260035|   1.554473|
#|-    |-        |                 0.70| 0.4409989| 0.2157869|      1.0767348|   1.826421|

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
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                              |Autres_(33%, 50%, 70%)                                                               |
#|:----|:----------------------------------------------------------------------------------|:------------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.00768847217963582, 0.00746025719887422, 0.00948874273516016 |Bias              , -0.0287619645011096 , -0.00247936584959296, 0.0194767289126385   |
#|-    |MCSE_var           , 0.00249312727815207, 0.00281292460561937, 0.00368024362570194 |Var               , 0.029556302228517   , 0.0278277187366774  , 0.0450181193470273   |
#|-    |MCSE_mse           , 0.00260251462196271, 0.00281195143507694, 0.00372133806848298 |Mse               , 0.0303244402260231  , 0.0277782105542201  , 0.0453074260774697   |
#|-    |%Cov               , 0.924              , 0.91               , 0.902               |Précision_var     , -176.071057961019   , -0.183104484824173  , -0.534199321462821   |

# nsim = 1000
#|n    |Erreur de Monte Carlo_(33%, 50%, 70%)                                              |Autres_(33%, 50%, 70%)                                                            |
#|:----|:----------------------------------------------------------------------------------|:---------------------------------------------------------------------------------|
#|1000 |MCSE_bias          , 0.0053329713236184 , 0.00515694737457127, 0.00682378115727317 |Bias               , -0.030906893936679 , -0.0131440549099147, 0.0203059332312411 |
#|-    |MCSE_var           , 0.0017783106257319 , 0.00194549783619299, 0.00264206055653868 |Var               , 0.0284405831385362, 0.0265941062240975, 0.0465639892823563    |
#|-    |MCSE_mse           , 0.00183975694599412, 0.00195018407251767, 0.00268874553674034 |Mse               , 0.0293673786482108, 0.0267402782973483, 0.0469297562174656    |
#|-    |%Cov , 0.915, 0.91 , 0.882                                                         |Précision_var     , -88.3605402613208 , -0.362640976372246, -0.872947850627278    |

