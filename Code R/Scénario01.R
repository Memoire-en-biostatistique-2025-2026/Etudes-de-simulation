
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

nsim <- 10 # Nombre de réplications

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
nsim <- 500
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
  resultats3[i,] <- TNDDR(dat, RN) # Estimateur doublement robuste
  
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

#   |n    |Methode |Erreur de Monte Carlo                      |Autres                                     |
#   |:----|:-------|:------------------------------------------|:------------------------------------------|
#   |1000 |RegLog  |MCSE_bias          , 0.00217543056861132   |Bias               , 0.00764669727604417   |
#   |-    |-       |MCSE_var            , 0.000206153565076863 |Var                 , 0.000206153565076863 |
#   |-    |-       |MCSE_mse            , 0.000216246580241676 |Mse                 , 0.000216246580241676 |
#   |-    |-       |%Cov , 0.943                               |Précision_var      , 0.00742922342628116   |

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

#   |n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
#   |:----|:-------|:------------------------------------------|:----------------------------------------|
#   |1000 |IPW     |MCSE_bias          , 0.00221880724759409   |Bias               , 0.00549099716653684 |
#   |-    |IPW     |MCSE_var            , 0.000211872879723978 |Var                , 0.00492310560197607 |
#   |-    |IPW     |MCSE_mse            , 0.000216246580241676 |Mse                , 0.00478623763992119 |
#   |-    |IPW     |%Cov , 0.921                               |Précision_var     , 0.0150996663466512   |

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
#|1000 |TNDDR   |MCSE_bias          , 0.0034136678808269    |Bias               , 0.000930345749421413  |
#|-    |-       |MCSE_var           , 0.000402232895985656  |Var                , 0.00582656420029462   |
#|-    |-       |MCSE_mse           , 0.000404130162849137  |Mse                , 0.00581577661510749   |
#|-    |-       |%Cov               , 0.932                 |Précision_var      , 0.00843627563555646   | # Avec CI 1
#|-    |-       |%Cov               , 0.932                 |Précision_var      , 0.0789107904397003    | # Avec CI 2
#|-    |-       |%Cov               , 0.932                 |Précision_var      , 0.0789107904397003    | # Avec CI 3

# Régression Lasso

#|n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
#|:----|:-------|:------------------------------------------|:----------------------------------------|
#|1000 |TNDDR   |MCSE_bias          , 0.00348818311854902   |Bias               , 0.0237237629440649  |
#|-    |-       |MCSE_var           , 0.000545147793575269  |Var                , 0.00608371073426517 |
#|-    |-       |MCSE_mse           , 0.000579060392633538  |Mse                , 0.00663436024102283 |
#|-    |-       |%Cov               , 0.74                  |Précision_var      , 0.222532880113912   | # Avec CI 1
#|-    |-       |%Cov               , 0.74                  |Précision_var      , 0.0651750510281717  | # Avec CI 2
#|-    |-       |%Cov               , 0.74                  |Précision_var      , 0.0651750510281717  | # Avec CI 3

# earth_GLM

#|n    |Methode |Erreur de Monte Carlo                    |Autres                                     |
#|:----|:-------|:----------------------------------------|:----------------------------------------  |
#|1000 |TNDDR   |MCSE_bias          , 0.00776555075993528 |Bias                 , -0.0264644191550205 |
#|-    |-       |MCSE_var           , 0.00292595961657473 |Var                  , 0.0301518893025657  |
#|-    |-       |MCSE_mse           , 0.0029290000237473  |Mse                  , 0.0307919510051732  |
#|-    |-       |%Cov               , 0.914               |Précision_var        , 0.677366581382007   | # Avec CI 1
#|-    |-       |%Cov               , 0.914               |Précision_var        , 0.677366581382007   | # Avec CI 2
#|-    |-       |%Cov               , 0.914               |Précision_var        , 0.677366581382007   | # Avec CI 3

# Réseaux de neurones

#|n    |Methode |Erreur de Monte Carlo                    |Autres                                    |
#|:----|:-------|:----------------------------------------|:--------------------------------------   |
#|1000 |TNDDR   |MCSE_bias          , 0.0133734403688143  Bias                  , -0.112985713207209 |
#|-    |-       |MCSE_var           , 0.00573928815306194 |Var                  , 0.0894244536491159 |
#|-    |-       |MCSE_mse           , 0.00478789261365004 |Mse                  , 0.102011376130759  |
#|-    |-       |%Cov               , 0.86                |Précision_var        , 0.586977884654825  | # Avec CI 1
#|-    |-       |%Cov               , 0.86                |Précision_var        , 0.586977884654825  | # Avec CI 2
#|-    |-       |%Cov               , 0.86                |Précision_var        , 0.586977884654825  | # Avec CI 3

## RN avec 1000 replications
#   |n    |Methode |Erreur de Monte Carlo                    |Autres                                   |
#   |:----|:-------|:----------------------------------------|:----------------------------------------|
#   |1000 |TNDDR   |MCSE_bias          , 0.00881153017434929 |Bias               , -0.0937663268196102 |
#   |-    |-       |MCSE_var           , 0.00373582186236799 |Var               , 0.0776430640134681   |
#   |-    |-       |MCSE_mse          , 0.0033289690698872   |Mse               , 0.0863575449946966   |
#   |-    |-       |%Cov , 0.861                             |Précision_var   , 1.80797415493798       |

# RN, 500 replications, plus de troncation
#   |n    |Methode |Erreur de Monte Carlo                    |Autres                                   |
#   |:----|:-------|:----------------------------------------|:----------------------------------------|
#   |1000 |TNDDR   |MCSE_bias          , 0.00777007509386557 |Bias               , -0.0170976621893679 |
#   |-    |-       |MCSE_var           , 0.00276525309482966 |Var              , 0.030187033482155     |
#   |-    |-       |MCSE_mse           , 0.00278017839460818 |Mse               , 0.0304189894675324   |
#   |-    |-       |%Cov , 0.914                             |Précision_var   , 1.15696399459045       |
