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
# Scénario 05 (co-infection ~ 30% dans l'échantillon): Sous-scénario 01 : W2_prev = 5%, 10% et 15% 
#                                                                         
# Pour W2_prev = 5% : co_inf_para1 = 1, co_inf_para2 = 3
# Pour W2_prev = 10% : 
# Pour W2_prev = 15% : 

################################################################################

################################################################################
################################################################################
# Pour W2_prev = 5% : co_inf_para1 = 1, co_inf_para2 = 3
# Proportion de co-infection ~30%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], W2_prev = 0.05,
                      co_inf_para1 = 1, co_inf_para2 = 3)
  
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
# 0.4764666
sd(l_vraiRRc)
# 0.005399389

l_vraiRRm

mean(l_vraiRRm)
# 0.4966849
sd(l_vraiRRm)
# 0.00536843

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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, W2_prev = 0.05,
                 co_inf_para1 = 1, co_inf_para2 = 3)
  
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
sum(rowSums(is.na(resultats2)) > 0)# 1

sum(rowSums(is.na(resultats3[,1:9])) > 0) #21
sum(rowSums(is.na(resultats3[,10:18])) > 0) #344
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
# |n    |Methode |Erreur de Monte Carlo                  |Autres                                 |
# |:----|:-------|:--------------------------------------|:--------------------------------------|
# |1000 |RegLog  |MCSE_bias         , 0.0119051556012247 |Bias             , 0.123255886325256   |
# |-    |-       |MCSE_var          , 0.0292567715038903 |Var              , 0.141732729889371   |
# |-    |-       |MCSE_mse          , 0.0314885161813705 |Mse              , 0.156783010673306   |
# |-    |-       |%Cov , 0.944                           |Précision_var     , 0.0226971533505271 |

################################ IPW ###########################################

Tab01$Methode <- c("IPW", "-", "-","-")

##    - statistiques descriptives

summary(resultats2$RRm)
mean(resultats2$RRm, na.rm = TRUE)
sd(resultats2$RRm, na.rm = TRUE)

mean(sqrt(na.omit(resultats2)$var_log_RRm))
sd(log(na.omit(resultats2$RRm)))

##    - biais, variance, moyenne de l'erreur-type, couverture des IC

resultats2$vrai_param <- rep(mean(l_vraiRRm), nsim) # vraie valeur du paramètre

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats2, RRm, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats2, RRm, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats2, RRm, vrai_param, criteria = "mse")

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
  list(name = "Précision_var", value = sd(log(na.omit(resultats2$RRm))) - mean(sqrt(na.omit(resultats2)$var_log_RRm)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode |Erreur de Monte Carlo                    |Autres                                   |
# |:----|:-------|:----------------------------------------|:----------------------------------------|
# |1000 |IPW     |MCSE_bias          , 0.00660970053837897 |Bias             , 0.034164864585566     |
# |-    |-       |MCSE_var           , 0.00152237812545456 |Var               , 0.0384892524034086   |
# |-    |-       |MCSE_mse           , 0.00170297128626108 |Mse               , 0.0396128022343517   |
# |-    |-       |%Cov            , 0.98410896708286       |Précision_var      , -0.0826328182408946 |
  
########################## TNDDR ###############################

Tab01$Methode <- c("TNDDR_RF", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_RF)
mean(resultats3$RRm_RF, na.rm = TRUE)
sd(resultats3$RRm_RF, na.rm =TRUE)

mean(sqrt(na.omit(resultats3)$`var_log_RRm-RF`))
sd(log(na.omit(resultats3)$RRm_RF))

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
  list(name = "Précision_var", value = sd(log(na.omit(resultats3)$RRm_RF)) - mean(sqrt(na.omit(resultats3)$`var_log_RRm-RF`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                   |
# |:----|:--------|:----------------------------------------|:----------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias         , 0.0100643712388581   |Bias               , -0.0802877062185647 |
# |-    |-        |MCSE_var           , 0.00332821075955121 |Var               , 0.0991644454964496   |
# |-    |-        |MCSE_mse           , 0.00320163385816734 |Mse              , 0.105509269697855     |
# |-    |-        |%Cov             , 0.979570990806946     |Précision_var    , 0.195575837197247     |
  
################################################################################

# Régression Lasso

Tab01$Methode <- c("TNDDR_Lasso", "-", "-","-")

##    - statistiques descriptives

summary(resultats3$RRm_Lasso)
mean(resultats3$RRm_Lasso, na.rm = TRUE)
sd(resultats3$RRm_Lasso, na.rm = TRUE)

mean(sqrt(na.omit(resultats3)$`var_log_RRm-Lasso`))
sd(log(na.omit(resultats3)$RRm_Lasso))

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
  list(name = "Précision_var", value = sd(log(na.omit(resultats3)$RRm_Lasso)) - mean(sqrt(na.omit(resultats3)$`var_log_RRm-Lasso`)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode     |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:-----------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias         , 0.0101796247762083   |Bias             , -0.09125779736688   |
# |-    |-           |MCSE_var           , 0.00340445109271167 |Var              , 0.067977842943363   |
# |-    |-           |MCSE_mse          , 0.0032413070855722   |Mse               , 0.0762022037630331 |
# |-    |-           |%Cov             , 0.844512195121951     |Précision_var    , 0.179090699086957   |
  
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
# |n    |Methode    |Erreur de Monte Carlo                    |Autres                                   |
# |:----|:----------|:----------------------------------------|:----------------------------------------|
# |1000 |TNDDR_Mars |MCSE_bias         , 0.0111022391846681   |Bias               , -0.0525983279926419 |
# |-    |-          |MCSE_var           , 0.00337789044064272 |Var             , 0.12325971491358       |
# |-    |-          |MCSE_mse           , 0.00328386292754602 |Mse              , 0.125903039306288     |
# |-    |-          |%Cov , 0.959                             |Précision_var    , -3.59759506742694     |

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
# |1000 |TNDDR_RN |MCSE_bias         , 0.0136288482575569   |Bias              , -0.155018335476325 |
# |-    |-        |MCSE_var           , 0.00476533767380971 |Var              , 0.185745504827512   |
# |-    |-        |MCSE_mse           , 0.00249874864452799 |Mse              , 0.209590443656535   |
# |-    |-        |%Cov , 0.871                             |Précision_var    , -4.39419207116346   |

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
# |1000 |TNDDR_GLM |MCSE_bias        , 0.010991860738328     |Bias              , -0.130844960853053 |
# |-    |-         |MCSE_var         , 0.003806227218257     |Var              , 0.120821002490796   |
# |-    |-         |MCSE_mse           , 0.00332897533130606 |Mse              , 0.137820585268942   |
# |-    |-         |%Cov , 0.927                             |Précision_var    , -11.1372549333478   |

################################################################################
################################################################################
# Pour I2_prev = 10% : co_inf_para1 = 3, co_inf_para2 = 4
# Proportion de co-infection ~30%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], I2_prev = 0.1,
                      co_inf_para1 = 3, co_inf_para2 = 4)
  
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
# 0.4074851
sd(l_vraiRRc)
# 0.01280726

l_vraiRRm

mean(l_vraiRRm)
# 0.4095705
sd(l_vraiRRm)
# 0.01280394

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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, I2_prev = 0.1,
                 co_inf_para1 = 3, co_inf_para2 = 4)
  
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
# |1000 |RegLog  |MCSE_bias          , 0.00260871026035455   |Bias             , 0.119368201729485     |
# |-    |-       |MCSE_var            , 0.000323631333741244 |Var                , 0.00680536922247909 |
# |-    |-       |MCSE_mse           , 0.00078087371080256   |Mse               , 0.0210473314373876   |
# |-    |-       |%Cov , 0.597                               |Précision_var     , 0.0105341459042426   |
  
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
# |1000 |IPW     |MCSE_bias          , 0.00259702207081865   |Bias             , 0.112001967943227     |
# |-    |-       |MCSE_var            , 0.000326574752520565 |Var                , 0.00674452363631919 |
# |-    |-       |MCSE_mse           , 0.00078087371080256   |Mse               , 0.0210473314373876   |
# |-    |-       |%Cov , 0.638                               |Précision_var      , 0.00991625633742094 |

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
# |1000 |TNDDR_RF |MCSE_bias          , 0.00269957031288337   |Bias             , 0.101522562461882     |
# |-    |-        |MCSE_var            , 0.000340624815775078 |Var                , 0.00728767987420123 |
# |-    |-        |MCSE_mse            , 0.000728873642817244 |Mse               , 0.0175872228831539   |
# |-    |-        |%Cov , 0.825                               |Précision_var     , 0.0072538422476357   |

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
# |1000 |TNDDR_Lasso |MCSE_bias         , 0.0026358178984253     |Bias             , 0.107422722898376     |
# |-    |-           |MCSE_var           , 0.00033849564122848   |Var                , 0.00694753599365914 |
# |-    |-           |MCSE_mse            , 0.000753200347393506 |Mse               , 0.0184802298525668   |
# |-    |-           |%Cov , 0.667                               |Précision_var     , 0.0306732232286172   |
  
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
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00482170057184637 |Bias              , 0.0937727019392964 |
# |-    |-          |MCSE_var           , 0.00170980180212321 |Var               , 0.0232487964045436 |
# |-    |-          |MCSE_mse           , 0.00180671574499136 |Mse               , 0.0320188672371352 |
# |-    |-          |%Cov , 0.892                             |Précision_var    , 0.722270403664815   |

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
# |1000 |TNDDR_RN |MCSE_bias          , 0.00598404688914692 |Bias              , 0.0696075204421452 |
# |-    |-        |MCSE_var           , 0.00218182676687493 |Var               , 0.0358088171715089 |
# |-    |-        |MCSE_mse           , 0.00210095498659638 |Mse               , 0.0406182152564411 |
# |-    |-        |%Cov , 0.887                             |Précision_var    , 0.638699220302498   |
  
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
# |n    |Methode   |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:---------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_GLM |MCSE_bias          , 0.00267380674641801   |Bias             , 0.109760926919426     |
# |-    |-         |MCSE_var            , 0.000355449501479047 |Var                , 0.00714924251719044 |
# |-    |-         |MCSE_mse            , 0.000782271897575258 |Mse               , 0.0191895543528848   |
# |-    |-         |%Cov , 0.786                               |Précision_var      , 0.00976804036018636 |

################################################################################
################################################################################
# Pour I2_prev = 20% : co_inf_para1 = 2, co_inf_para2 = 3
# Proportion de co-infection ~40%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], I2_prev = 0.2,
                      co_inf_para1 = 2, co_inf_para2 = 3)
  
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
# 0.4166192
sd(l_vraiRRc)
# 0.01037789

l_vraiRRm

mean(l_vraiRRm)
# 0.42059
sd(l_vraiRRm)
# 0.01037681

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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, I2_prev = 0.2,
                 co_inf_para1 = 2, co_inf_para2 = 3)
  
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
# |1000 |RegLog  |MCSE_bias          , 0.00282994104418464   |Bias             , 0.121601530415067     |
# |-    |-       |MCSE_var            , 0.000370259871057647 |Var                , 0.00800856631356086 |
# |-    |-       |MCSE_mse            , 0.000866942372374566 |Mse               , 0.0227874899465338   |
# |-    |-       |%Cov , 0.643                               |Précision_var     , 0.0101923592583766   |

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
# |1000 |IPW     |MCSE_bias         , 0.0028011189453454     |Bias            , 0.10555342139345       |
# |-    |-       |MCSE_var            , 0.000359851976441348 |Var                , 0.00784626734597292 |
# |-    |-       |MCSE_mse            , 0.000866942372374566 |Mse               , 0.0227874899465338   |
# |-    |-       |%Cov , 0.728                               |Précision_var      , 0.00802647012432042 |

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
# |1000 |TNDDR_RF |MCSE_bias          , 0.00299421497408163   |Bias             , 0.087610451048832     |
# |-    |-        |MCSE_var            , 0.000470394371996655 |Var                , 0.00896532331101465 |
# |-    |-        |MCSE_mse            , 0.000749784726444859 |Mse               , 0.0166319491206834   |
# |-    |-        |%Cov , 0.894                               |Précision_var     , 0.0152630049390839   |

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
# |1000 |TNDDR_Lasso |MCSE_bias         , 0.0028175453875194     |Bias              , 0.0951678142060868   |
# |-    |-           |MCSE_var            , 0.000357217790252782 |Var                , 0.00793856201073183 |
# |-    |-           |MCSE_mse            , 0.000717022952409307 |Mse               , 0.0169875363094854   |
# |-    |-           |%Cov , 0.717                               |Précision_var     , 0.0408618224496662   |

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
# |1000 |TNDDR_Mars |MCSE_bias        , 0.005488260126024     |Bias              , 0.0672926195995585 |
# |-    |-          |MCSE_var           , 0.00194112100540389 |Var              , 0.030120999210905   |
# |-    |-          |MCSE_mse           , 0.00202357331714353 |Mse              , 0.034619174864265   |
# |-    |-          |%Cov, 0.95                               |Précision_var    , 0.587319712524824   |

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
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                     |
# |:----|:--------|:----------------------------------------|:------------------------------------------|
# |1000 |TNDDR_RN |MCSE_bias          , 0.00842685191140553 |Bias                , -0.00304690134784674 |
# |-    |-        |MCSE_var           , 0.00277958320089358 |Var               , 0.0710118331367591     |
# |-    |-        |MCSE_mse           , 0.00278255969669293 |Mse               , 0.0709501049114458     |
# |-    |-        |%Cov , 0.937                             |Précision_var   , 1.33302482248194         |

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
# |n    |Methode   |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:---------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_GLM |MCSE_bias          , 0.00283935693601321   |Bias             , 0.096449748128369     |
# |-    |-         |MCSE_var            , 0.000393766601426405 |Var                , 0.00806194781008632 |
# |-    |-         |MCSE_mse            , 0.000756793202283021 |Mse              , 0.017356439776302     |
# |-    |-         |%Cov , 0.886                               |Précision_var      , 0.00202892000316651 |
