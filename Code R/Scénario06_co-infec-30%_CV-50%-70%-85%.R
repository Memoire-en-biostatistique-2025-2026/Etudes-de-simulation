
# Execute R code directly from the following files

source("Code R/Génération des données.R")
source("Code R/Fonc01_RegLog.R")
source("Code R/Fonc02_IPW.R")
source("Code R/Fonc03_TNDDR.R")

# Install and load required packages

library(simhelpers)
library(dplyr)
library(tibble)
library(knitr)
library(dplyr)
library(kableExtra)

################################################################################
# Scénario 06 (co-infection ~ 30% dans l'échantillon): Sous-scénario 01 : Couverture vaccinale = 50%, 70% et 85% 
#                                                                         I1_prev = 15% et 
#                                                                         I2 = 50%
# Pour CV = 50% : co_inf_para1 = 0.5, co_inf_para2 = 0.5
# Pour CV = 70% : co_inf_para1 = 0.5, co_inf_para2 = 0.5
# Pour CV = 85% : co_inf_para1 = 0.5, co_inf_para2 = 0.5

################################################################################

################################################################################
################################################################################
# Pour CV = 50% : co_inf_para1 = 0.5, co_inf_para2 = 0.5
# Proportion de co-infection ~30%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], CV = 0.5,
                      co_inf_para1 = 0.5, co_inf_para2 = 0.5)
  
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
# 0.4126738
sd(l_vraiRRc)
# 0.004989355

l_vraiRRm

mean(l_vraiRRm)
# 0.4237495
sd(l_vraiRRm)
# 0.005026357

################################################################################
# Initialize objects to store results

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

################################################################################
########################## Analysis of the results #############################

nsim <- 1000 # Number of replicas

# Initialize objects to store results

resultats <- data.frame(matrix(ncol = 4, 
                               nrow = nsim))

colnames(resultats) <- c("coe_reg", # Logistic regression coefficient
                         "err_reg", # Standard error of the coefficient
                         "RRc", # Conditional risk ratio
                         "est_VE") # ^VE (vaccine effectiveness)

resultats2 <- data.frame(matrix(ncol = 5, 
                                nrow = nsim))

colnames(resultats2) <- c("RRm",# Marginal relative risk
                          "VE",# ^VE (vaccine effectiveness)
                          "var_log_RRm",# Variance of the log of mRR
                          "IC_inf", # Lower bound of the confidence interval
                          "IC_sup") # Upper bound of the confidence interval

resultats3 <- data.frame(matrix(ncol = 45, 
                                nrow = nsim))

colnames(resultats3) <- c("RRm_RF", "VE_RF", "var_log_RRm-RF", "IC_inf1-RF", "IC_sup1-RF", "IC_inf2-RF", "IC_sup2-RF", 
                          "IC_inf3-RF", "IC_sup3-RF",# Random Forest
                          
                          "RRm_Lasso", "VE_Lasso","var_log_RRm-Lasso", "IC_inf1-Lasso", "IC_sup1-Lasso", 
                          "IC_inf2-Lasso", "IC_sup2-Lasso", 
                          "IC_inf3-Lasso", "IC_sup3-Lasso",# Lasso regression
                          
                          "RRm_Mars", "VE_Mars","var_log_RRm-Mars", "IC_inf1-Mars", "IC_sup1-Mars", 
                          "IC_inf2-Mars", "IC_sup2-Mars", 
                          "IC_inf3-Mars", "IC_sup3-Mars",# MARS
                          
                          "RRm_RN", "VE_RN","var_log_RRm-RN", "IC_inf1-RN", "IC_sup1-RN", 
                          "IC_inf2-RN", "IC_sup2-RN", 
                          "IC_inf3-RN", "IC_sup3-RN",# Neural networks
                          
                          "RRm_GLM", "VE_GLM","var_log_RRm-GLM", "IC_inf1-GLM", "IC_sup1-GLM", 
                          "IC_inf2-GLM", "IC_sup2-GLM", 
                          "IC_inf3-GLM", "IC_sup3-GLM"# GLM
                          
                          
                          # Marginal relative risk  
                          # ^VE (vaccine effectiveness)
                          # Variance of the log of the mRR
                          # Lower bound of the first confidence interval
                          # Its upper bound
                          # Lower bound of the second confidence interval
                          # Its upper bound
                          # Lower bound of the third confidence interval
                          # Its upper bound
                          
)

methode <- list(RandomForest, Lasso, Mars, RN, PM)

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, CV = 0.5,
                 co_inf_para1 = 0.5, co_inf_para2 = 0.5)
  
  tryCatch({
    resultats[i,] <- RegLog(dat) # Logistic regression
  }, error = function(e){})
  
  tryCatch({
    resultats2[i,] <- IPW(dat)   # IPW
  }, error = function(e){})
  
  l <- list()
  
  for(j in methode) {
    k <- rep(NA, 9)
    tryCatch({
      
      k <- TNDDR(dat, j) # List of results for method j
      
    }, error = function(e){})
    
    l <- append(l, k) # Combining the results of different methods for TNDDR
    
    resultats3[i,] <- l
    
    
  }
  
  # Add a line to track progress
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}

# How many replicates with NA?

sum(rowSums(is.na(resultats)) > 0)
sum(rowSums(is.na(resultats2)) > 0)

sum(rowSums(is.na(resultats3[,1:9])) > 0)
sum(rowSums(is.na(resultats3[,10:18])) > 0)
sum(rowSums(is.na(resultats3[,19:27])) > 0)
sum(rowSums(is.na(resultats3[,28:36])) > 0)
sum(rowSums(is.na(resultats3[,37:45])) > 0)

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
# |1000 |RegLog  |MCSE_bias          , 0.00328442525317894   |Bias              , 0.0810017629345915   |
# |-    |-       |MCSE_var            , 0.000556080134958981 |Var               , 0.0107874492437195   |
# |-    |-       |MCSE_mse           , 0.00089199939527236   |Mse               , 0.0173379473929876   |
# |-    |-       |%Cov , 0.876                               |Précision_var      , 0.00720292746937184 |

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
# |1000 |IPW     |MCSE_bias         , 0.0034951313358671     |Bias              , 0.0538988892837293   |
# |-    |-       |MCSE_var            , 0.000576573634202178 |Var               , 0.0122159430549602   |
# |-    |-       |MCSE_mse           , 0.00089199939527236   |Mse               , 0.0173379473929876   |
# |-    |-       |%Cov , 0.923                               |Précision_var      , 0.00793962481856178 |

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
# |1000 |TNDDR_RF |MCSE_bias          , 0.00384037571760654   |Bias              , 0.0671380201355453 |
# |-    |-        |MCSE_var            , 0.000890282212025785 |Var              , 0.014748485652382   |
# |-    |-        |MCSE_mse           , 0.00117675430450768   |Mse               , 0.0192412509144505 |
# |-    |-        |%Cov , 0.969                               |Précision_var    , 0.109338080327853   |

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
# |n    |Methode     |Erreur de Monte Carlo                      |Autres                                 |
# |:----|:-----------|:------------------------------------------|:--------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias          , 0.00350538845619588   |Bias             , 0.044866823742757   |
# |-    |-           |MCSE_var            , 0.000593721373047661 |Var               , 0.0122877482288313 |
# |-    |-           |MCSE_mse            , 0.000749348020728449 |Mse               , 0.0142884923533661 |
# |-    |-           |%Cov, 0.72                                 |Précision_var    , 0.110293882533244   |

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
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00693601797593869 |Bias              , 0.0295498976121155 |
# |-    |-          |MCSE_var           , 0.00238312954301461 |Var               , 0.0481083453625446 |
# |-    |-          |MCSE_mse          , 0.0024077203677201   |Mse               , 0.0489334334660686 |
# |-    |-          |%Cov , 0.934                             |Précision_var    , 0.754431587411811   |

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
# |1000 |TNDDR_RN |MCSE_bias         , 0.0106946149427441   |Bias              , -0.108883227465173 |
# |-    |-        |MCSE_var           , 0.00406813015528204 |Var              , 0.114374788773565   |
# |-    |-        |MCSE_mse           , 0.00313462725473016 |Mse              , 0.126115971208024   |
# |-    |-        |%Cov , 0.774                             |Précision_var   , -17.973800644942     |

################################################################################

# GLM

Tab01$Methode <- c("TNDDR_GLM", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_GLM)
mean(resultats3$RRm_GLM)
sd(resultats3$RRm_GLM)

mean(sqrt(resultats3$`var_log_RRm-GLM`))
sd(log(resultats3$RRm_GLM))

help("calc_absolute") # calculate various performance metrics

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_GLM, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_GLM, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_GLM, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

coverage <- calc_coverage(resultats3, `IC_inf2-GLM`, `IC_sup2-GLM`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_GLM)) - mean(sqrt(resultats3$`var_log_RRm-GLM`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode   |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:---------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_GLM |MCSE_bias          , 0.00431971371882887 |Bias              , 0.0692894654627735 |
# |-    |-         |MCSE_var           , 0.00116403095987496 |Var               , 0.0186599266126384 |
# |-    |-         |MCSE_mse        , 0.00152961202535       |Mse               , 0.0234422967101426 |
# |-    |-         |%Cov , 0.964                             |Précision_var     , 0.0956901852445269 |

################################################################################
################################################################################
# Pour CV = 70% : co_inf_para1 = 0.5, co_inf_para2 = 0.5
# Proportion de co-infection ~30%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], CV = 0.7,
                      co_inf_para1 = 0.5, co_inf_para2 = 0.5)
  
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
# 0.4243789
sd(l_vraiRRc)
# 0.005631539

l_vraiRRm

mean(l_vraiRRm)
# 0.4372749
sd(l_vraiRRm)
# 0.00566205

################################################################################
# Initialize objects to store results

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

################################################################################
########################## Analysis of the results #############################

nsim <- 1000 # Number of replicas

# Initialize objects to store results

resultats <- data.frame(matrix(ncol = 4, 
                               nrow = nsim))

colnames(resultats) <- c("coe_reg", # Logistic regression coefficient
                         "err_reg", # Standard error of the coefficient
                         "RRc", # Conditional risk ratio
                         "est_VE") # ^VE (vaccine effectiveness)

resultats2 <- data.frame(matrix(ncol = 5, 
                                nrow = nsim))

colnames(resultats2) <- c("RRm",# Marginal relative risk
                          "VE",# ^VE (vaccine effectiveness)
                          "var_log_RRm",# Variance of the log of mRR
                          "IC_inf", # Lower bound of the confidence interval
                          "IC_sup") # Upper bound of the confidence interval

resultats3 <- data.frame(matrix(ncol = 45, 
                                nrow = nsim))

colnames(resultats3) <- c("RRm_RF", "VE_RF", "var_log_RRm-RF", "IC_inf1-RF", "IC_sup1-RF", "IC_inf2-RF", "IC_sup2-RF", 
                          "IC_inf3-RF", "IC_sup3-RF",# Random Forest
                          
                          "RRm_Lasso", "VE_Lasso","var_log_RRm-Lasso", "IC_inf1-Lasso", "IC_sup1-Lasso", 
                          "IC_inf2-Lasso", "IC_sup2-Lasso", 
                          "IC_inf3-Lasso", "IC_sup3-Lasso",# Lasso regression
                          
                          "RRm_Mars", "VE_Mars","var_log_RRm-Mars", "IC_inf1-Mars", "IC_sup1-Mars", 
                          "IC_inf2-Mars", "IC_sup2-Mars", 
                          "IC_inf3-Mars", "IC_sup3-Mars",# MARS
                          
                          "RRm_RN", "VE_RN","var_log_RRm-RN", "IC_inf1-RN", "IC_sup1-RN", 
                          "IC_inf2-RN", "IC_sup2-RN", 
                          "IC_inf3-RN", "IC_sup3-RN",# Neural networks
                          
                          "RRm_GLM", "VE_GLM","var_log_RRm-GLM", "IC_inf1-GLM", "IC_sup1-GLM", 
                          "IC_inf2-GLM", "IC_sup2-GLM", 
                          "IC_inf3-GLM", "IC_sup3-GLM"# GLM
                          
                          
                          # Marginal relative risk  
                          # ^VE (vaccine effectiveness)
                          # Variance of the log of the mRR
                          # Lower bound of the first confidence interval
                          # Its upper bound
                          # Lower bound of the second confidence interval
                          # Its upper bound
                          # Lower bound of the third confidence interval
                          # Its upper bound
                          
)

methode <- list(RandomForest, Lasso, Mars, RN, PM)

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, CV = 0.7,
                 co_inf_para1 = 0.5, co_inf_para2 = 0.5)
  
  tryCatch({
    resultats[i,] <- RegLog(dat) # Logistic regression
  }, error = function(e){})
  
  tryCatch({
    resultats2[i,] <- IPW(dat)   # IPW
  }, error = function(e){})
  
  l <- list()
  
  for(j in methode) {
    k <- rep(NA, 9)
    tryCatch({
      
      k <- TNDDR(dat, j) # List of results for method j
      
    }, error = function(e){})
    
    l <- append(l, k) # Combining the results of different methods for TNDDR
    
    resultats3[i,] <- l
    
    
  }
  
  # Add a line to track progress
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}

# How many replicates with NA?

sum(rowSums(is.na(resultats)) > 0)
sum(rowSums(is.na(resultats2)) > 0)

sum(rowSums(is.na(resultats3[,1:9])) > 0)
sum(rowSums(is.na(resultats3[,10:18])) > 0) # 2
sum(rowSums(is.na(resultats3[,19:27])) > 0)
sum(rowSums(is.na(resultats3[,28:36])) > 0)
sum(rowSums(is.na(resultats3[,37:45])) > 0)

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
# |1000 |RegLog  |MCSE_bias          , 0.00354123775557352   |Bias              , 0.0627651391129958   |
# |-    |-       |MCSE_var            , 0.000610671256945237 |Var               , 0.0125403648414994   |
# |-    |-       |MCSE_mse            , 0.000878095874259401 |Mse               , 0.0164672871645316   |
# |-    |-       |%Cov , 0.901                               |Précision_var      , 0.00264523738571165 |

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
# |1000 |IPW     |MCSE_bias          , 0.00386210218152429   |Bias              , 0.0257236505598142   |
# |-    |-       |MCSE_var            , 0.000828343071096757 |Var               , 0.0149158332605347   |
# |-    |-       |MCSE_mse            , 0.000878095874259401 |Mse               , 0.0164672871645316   |
# |-    |-       |%Cov , 0.949                               |Précision_var      , 0.00330176712403529 |

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
# |1000 |TNDDR_RF |MCSE_bias          , 0.00562343388872076 |Bias              , 0.0884907129141091 |
# |-    |-        |MCSE_var           , 0.00182444198128456 |Var               , 0.0316230087008131 |
# |-    |-        |MCSE_mse          , 0.0022487606256927   |Mse               , 0.0394219919641596 |
# |-    |-        |%Cov , 0.965                             |Précision_var     , -0.382711620383291 |

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
# |n    |Methode     |Erreur de Monte Carlo                      |Autres                                 |
# |:----|:-----------|:------------------------------------------|:--------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias          , 0.00391049315039648   |Bias              , 0.0209138420001072 |
# |-    |-           |MCSE_var            , 0.000762395286547975 |Var               , 0.0152919566792978 |
# |-    |-           |MCSE_mse            , 0.000817396619136088 |Mse              , 0.015714053509824   |
# |-    |-           |%Cov , 0.549                               |Précision_var    , 0.186008081717384   |

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
# |1000 |TNDDR_Mars |MCSE_bias         , 0.0076050877592958   |Bias              , 0.0507910567674996 |
# |-    |-          |MCSE_var          , 0.0026516981955901   |Var               , 0.0578373598265908 |
# |-    |-          |MCSE_mse           , 0.00273766193767428 |Mse               , 0.0603592539143235 |
# |-    |-          |%Cov , 0.895                             |Précision_var    , 0.837180643029826   |

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
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                               |
# |:----|:--------|:----------------------------------------|:------------------------------------|
# |1000 |TNDDR_RN |MCSE_bias         , 0.0109004340868135   |Bias             , -0.18337572638961 |
# |-    |-        |MCSE_var           , 0.00501021833006288 |Var              , 0.118819463280966 |
# |-    |-        |MCSE_mse           , 0.00289297606246073 |Mse              , 0.152327300846603 |
# |-    |-        |%Cov , 0.557                             |Précision_var    , -1.14206207233916 |

################################################################################

# GLM

Tab01$Methode <- c("TNDDR_GLM", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_GLM)
mean(resultats3$RRm_GLM)
sd(resultats3$RRm_GLM)

mean(sqrt(resultats3$`var_log_RRm-GLM`))
sd(log(resultats3$RRm_GLM))

help("calc_absolute") # calculate various performance metrics

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_GLM, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_GLM, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_GLM, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

coverage <- calc_coverage(resultats3, `IC_inf2-GLM`, `IC_sup2-GLM`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_GLM)) - mean(sqrt(resultats3$`var_log_RRm-GLM`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode   |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:---------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_GLM |MCSE_bias          , 0.00617728791393877 |Bias             , 0.082891781968422   |
# |-    |-         |MCSE_var           , 0.00206103473257261 |Var              , 0.038158885971694   |
# |-    |-         |MCSE_mse           , 0.00249697410777121 |Mse               , 0.0449917746036227 |
# |-    |-         |%Cov , 0.951                             |Précision_var     , 0.0540649457740704 |

################################################################################
################################################################################
# Pour CV = 85% : co_inf_para1 = 0.5, co_inf_para2 = 0.5
# Proportion de co-infection ~30%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], CV = 0.85,
                      co_inf_para1 = 0.5, co_inf_para2 = 0.5)
  
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
# 0.4285399
sd(l_vraiRRc)
# 0.003754405

l_vraiRRm

mean(l_vraiRRm)
# 0.4430931
sd(l_vraiRRm)
# 0.00373983

################################################################################
# Initialize objects to store results

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

################################################################################
########################## Analysis of the results #############################

nsim <- 1000 # Number of replicas

# Initialize objects to store results

resultats <- data.frame(matrix(ncol = 4, 
                               nrow = nsim))

colnames(resultats) <- c("coe_reg", # Logistic regression coefficient
                         "err_reg", # Standard error of the coefficient
                         "RRc", # Conditional risk ratio
                         "est_VE") # ^VE (vaccine effectiveness)

resultats2 <- data.frame(matrix(ncol = 5, 
                                nrow = nsim))

colnames(resultats2) <- c("RRm",# Marginal relative risk
                          "VE",# ^VE (vaccine effectiveness)
                          "var_log_RRm",# Variance of the log of mRR
                          "IC_inf", # Lower bound of the confidence interval
                          "IC_sup") # Upper bound of the confidence interval

resultats3 <- data.frame(matrix(ncol = 45, 
                                nrow = nsim))

colnames(resultats3) <- c("RRm_RF", "VE_RF", "var_log_RRm-RF", "IC_inf1-RF", "IC_sup1-RF", "IC_inf2-RF", "IC_sup2-RF", 
                          "IC_inf3-RF", "IC_sup3-RF",# Random Forest
                          
                          "RRm_Lasso", "VE_Lasso","var_log_RRm-Lasso", "IC_inf1-Lasso", "IC_sup1-Lasso", 
                          "IC_inf2-Lasso", "IC_sup2-Lasso", 
                          "IC_inf3-Lasso", "IC_sup3-Lasso",# Lasso regression
                          
                          "RRm_Mars", "VE_Mars","var_log_RRm-Mars", "IC_inf1-Mars", "IC_sup1-Mars", 
                          "IC_inf2-Mars", "IC_sup2-Mars", 
                          "IC_inf3-Mars", "IC_sup3-Mars",# MARS
                          
                          "RRm_RN", "VE_RN","var_log_RRm-RN", "IC_inf1-RN", "IC_sup1-RN", 
                          "IC_inf2-RN", "IC_sup2-RN", 
                          "IC_inf3-RN", "IC_sup3-RN",# Neural networks
                          
                          "RRm_GLM", "VE_GLM","var_log_RRm-GLM", "IC_inf1-GLM", "IC_sup1-GLM", 
                          "IC_inf2-GLM", "IC_sup2-GLM", 
                          "IC_inf3-GLM", "IC_sup3-GLM"# GLM
                          
                          
                          # Marginal relative risk  
                          # ^VE (vaccine effectiveness)
                          # Variance of the log of the mRR
                          # Lower bound of the first confidence interval
                          # Its upper bound
                          # Lower bound of the second confidence interval
                          # Its upper bound
                          # Lower bound of the third confidence interval
                          # Its upper bound
                          
)

methode <- list(RandomForest, Lasso, Mars, RN, PM)

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, CV = 0.85,
                 co_inf_para1 = 0.5, co_inf_para2 = 0.5)
  
  tryCatch({
    resultats[i,] <- RegLog(dat) # Logistic regression
  }, error = function(e){})
  
  tryCatch({
    resultats2[i,] <- IPW(dat)   # IPW
  }, error = function(e){})
  
  l <- list()
  
  for(j in methode) {
    k <- rep(NA, 9)
    tryCatch({
      
      k <- TNDDR(dat, j) # List of results for method j
      
    }, error = function(e){})
    
    l <- append(l, k) # Combining the results of different methods for TNDDR
    
    resultats3[i,] <- l
    
    
  }
  
  # Add a line to track progress
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
}

# How many replicates with NA?

sum(rowSums(is.na(resultats)) > 0)
sum(rowSums(is.na(resultats2)) > 0)

sum(rowSums(is.na(resultats3[,1:9])) > 0) 
sum(rowSums(is.na(resultats3[,10:18])) > 0) # 68
sum(rowSums(is.na(resultats3[,19:27])) > 0)
sum(rowSums(is.na(resultats3[,28:36])) > 0)
sum(rowSums(is.na(resultats3[,37:45])) > 0)

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
# |n    |Methode |Erreur de Monte Carlo                    |Autres                                   |
# |:----|:-------|:----------------------------------------|:----------------------------------------|
# |1000 |RegLog  |MCSE_bias          , 0.00469476834981632 |Bias              , 0.0654172587599151   |
# |-    |-       |MCSE_var           , 0.00122873936545015 |Var              , 0.022040849858437     |
# |-    |-       |MCSE_mse           , 0.00157961420276785 |Mse               , 0.0262982267522403   |
# |-    |-       |%Cov , 0.918                             |Précision_var      , 0.00890590357553944 |

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
# |n    |Methode |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:-------|:----------------------------------------|:--------------------------------------|
# |1000 |IPW     |MCSE_bias          , 0.00523749562850913 |Bias              , 0.0271377813327935 |
# |-    |-       |MCSE_var           , 0.00147288090711237 |Var               , 0.0274313604586522 |
# |-    |-       |MCSE_mse           , 0.00157961420276785 |Mse               , 0.0262982267522403 |
# |-    |-       |%Cov , 0.948                             |Précision_var     , 0.0155419273179212 |

########################## TNDDR ###############################

Tab01$Methode <- c("TNDDR_RF", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_RF)
mean(resultats3$RRm_RF, na.rm = TRUE)
sd(resultats3$RRm_RF, na.rm = TRUE)

mean(sqrt(na.omit(resultats3$`var_log_RRm-RF`)))
sd(log(na.omit(resultats3$RRm_RF)))

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
  list(name = "Précision_var", value = sd(log(na.omit(resultats3$RRm_RF))) - mean(sqrt(na.omit(resultats3$`var_log_RRm-RF`))) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:--------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00962076269589796 |Bias              , 0.0923212713849915 |
# |-    |-        |MCSE_var          , 0.0032363330268851   |Var               , 0.0925590748507817 |
# |-    |-        |MCSE_mse          , 0.0035048973550282   |Mse              , 0.100989732926072   |
# |-    |-        |%Cov , 0.838                             |Précision_var    , 0.792849798123044   |

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
# |n    |Methode     |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:-----------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias          , 0.00532134577317481 |Bias             , 0.027099096362143   |
# |-    |-           |MCSE_var           , 0.00134164956489766 |Var               , 0.0263911838207228 |
# |-    |-           |MCSE_mse           , 0.00147113599474243 |Mse               , 0.0270972281235299 |
# |-    |-           |%Cov             , 0.452789699570815     |Précision_var    , 0.265943275361501   |

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
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00901101369034116 |Bias              , 0.0342708918731162 |
# |-    |-          |MCSE_var           , 0.00302321482648016 |Var               , 0.0811983677275158 |
# |-    |-          |MCSE_mse           , 0.00306771666700249 |Mse               , 0.0822916633895671 |
# |-    |-          |%Cov , 0.799                             |Précision_var    , 0.501357288732908   |

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
# |1000 |TNDDR_RN |MCSE_bias         , 0.0105273638400162   |Bias              , -0.239288402106295 |
# |-    |-        |MCSE_var           , 0.00576994421850924 |Var              , 0.110825389420081   |
# |-    |-        |MCSE_mse           , 0.00254456133378121 |Mse              , 0.167973503413245   |
# |-    |-        |%Cov , 0.387                             |Précision_var    , -3.45325362000664   |

################################################################################

# GLM

Tab01$Methode <- c("TNDDR_GLM", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_GLM)
mean(resultats3$RRm_GLM)
sd(resultats3$RRm_GLM)

mean(sqrt(resultats3$`var_log_RRm-GLM`))
sd(log(resultats3$RRm_GLM))

help("calc_absolute") # calculate various performance metrics

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_GLM, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_GLM, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_GLM, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

coverage <- calc_coverage(resultats3, `IC_inf2-GLM`, `IC_sup2-GLM`, vrai_param)

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_GLM)) - mean(sqrt(resultats3$`var_log_RRm-GLM`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode   |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:---------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_GLM |MCSE_bias          , 0.00967984991417517 |Bias              , 0.0655980632361706 |
# |-    |-         |MCSE_var           , 0.00323836298001078 |Var               , 0.0936994943609569 |
# |-    |-         |MCSE_mse         , 0.003449300650375     |Mse               , 0.0979089007669326 |
# |-    |-         |%Cov , 0.809                             |Précision_var    , 0.681121811511624   |