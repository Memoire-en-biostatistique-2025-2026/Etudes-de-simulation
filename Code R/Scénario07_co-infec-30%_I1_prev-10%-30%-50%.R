
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
# Scénario 07 (co-infection ~ 30% dans l'échantillon): Sous-scénario 01 : Couverture vaccinale = 33% 
#                                                                         I1_prev = 10%, 30% et 50% et 
#                                                                         I2 = 50%
# Pour I1_rev = 10% : co_inf_para1 = 1, co_inf_para2 = 2
# Pour I1_rev = 30% : co_inf_para1 = 1, co_inf_para2 = -1
# Pour I1_rev = 50% : co_inf_para1 = 1.5, co_inf_para2 = -1.5

################################################################################

################################################################################
################################################################################
# Pour I1_prev = 10% : co_inf_para1 = 1, co_inf_para2 = 2
# Proportion de co-infection ~30%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], I1_prev = 0.10,
                      co_inf_para1 = 1, co_inf_para2 = 2)
  
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
# 0.4239555
sd(l_vraiRRc)
# 0.005397573

l_vraiRRm

mean(l_vraiRRm)
# 0.4330938
sd(l_vraiRRm)
# 0.005406291

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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, I1_prev = 0.1,
                 co_inf_para1 = 1, co_inf_para2 = 2)
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
sum(rowSums(is.na(resultats3[,10:18])) > 0) # 44
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
# |n    |Methode |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:-------|:----------------------------------------|:--------------------------------------|
# |1000 |RegLog  |MCSE_bias          , 0.00628309116985792 |Bias             , 0.101727803891357   |
# |-    |-       |MCSE_var           , 0.00356792444605976 |Var               , 0.0394772346487466 |
# |-    |-       |MCSE_mse           , 0.00442717519793174 |Mse               , 0.0497863034986563 |
# |-    |-       |%Cov, 0.93                               |Précision_var     , 0.0193522785525343 |
  
################################ IPW ###########################################

Tab01$Methode <- c("IPW", "-", "-","-")

##    - statistiques descriptives

summary(resultats2$RRm)
mean(resultats2$RRm, na.rm = TRUE)
sd(resultats2$RRm, na.rm = TRUE)

mean(sqrt(na.omit(resultats2)$var_log_RRm))
sd(log(na.omit(resultats2)$RRm))

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
  list(name = "Précision_var", value = sd(log(na.omit(resultats2)$RRm)) - mean(sqrt(na.omit(resultats2)$var_log_RRm)) # Voir si la variance est bien estimée
       
  )
  
)

kable(Tab01)
# |n    |Methode |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:-------|:----------------------------------------|:--------------------------------------|
# |1000 |IPW     |MCSE_bias         , 0.0090739123542567   |Bias             , 0.108774276222893   |
# |-    |-       |MCSE_var          , 0.0216393234184862   |Var               , 0.0823358854127323 |
# |-    |-       |MCSE_mse           , 0.00442717519793174 |Mse               , 0.0497863034986563 |
# |-    |-       |%Cov , 0.942                             |Précision_var     , 0.0309303314511636 |

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
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                     |
# |:----|:--------|:----------------------------------------|:------------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00794106404482242 |Bias                , -0.00678488115211556 |
# |-    |-        |MCSE_var           , 0.00269225360927307 |Var               , 0.0630604981639713     |
# |-    |-        |MCSE_mse           , 0.00267093369759151 |Mse               , 0.0630434722780557     |
# |-    |-        |%Cov , 0.973                             |Précision_var     , -0.467228230030712     |
  
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
# |1000 |TNDDR_Lasso |MCSE_bias         , 0.0104349805031381   |Bias              , -0.132341327349331 |
# |-    |-           |MCSE_var           , 0.00380741927237039 |Var              , 0.104097710104434   |
# |-    |-           |MCSE_mse           , 0.00291175193613568 |Mse              , 0.121503048210916   |
# |-    |-           |%Cov             , 0.993723849372385     |Précision_var    , -1.75271221737612   |

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
# |n    |Methode    |Erreur de Monte Carlo                    |Autres                                     |
# |:----|:----------|:----------------------------------------|:------------------------------------------|
# |1000 |TNDDR_Mars |MCSE_bias        , 0.010126818865004     |Bias                , -0.00554603876917303 |
# |-    |-          |MCSE_var           , 0.00325907048849278 |Var              , 0.102552460324601       |
# |-    |-          |MCSE_mse           , 0.00322970432871911 |Mse              , 0.102480666410305       |
# |-    |-          |%Cov , 0.953                             |Précision_var    , -5.80488013777676       |

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
# |1000 |TNDDR_RN |MCSE_bias         , 0.0128114036525453   |Bias              , -0.145303980054038 |
# |-    |-        |MCSE_var           , 0.00556684391921763 |Var              , 0.164132063548452   |
# |-    |-        |MCSE_mse           , 0.00277244569017329 |Mse              , 0.185081178104448   |
# |-    |-        |%Cov , 0.817                             |Précision_var    , -1.25753878021272   |

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
# |1000 |TNDDR_GLM |MCSE_bias         , 0.0108812718445013   |Bias              , -0.155005045470538 |
# |-    |-         |MCSE_var          , 0.0046135507818836   |Var              , 0.118402076953937   |
# |-    |-         |MCSE_mse           , 0.00294027932358262 |Mse              , 0.142310238998307   |
# |-    |-         |%Cov , 0.999                             |Précision_var    , -5.59110180839257   |

################################################################################
################################################################################
# Pour I1_prev = 30% : co_inf_para1 = 1, co_inf_para2 = -1
# Proportion de co-infection ~30%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], I1_prev = 0.3,
                      co_inf_para1 = 1, co_inf_para2 = -1)
  
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
# 0.4144212
sd(l_vraiRRc)
# 0.007771347

l_vraiRRm

mean(l_vraiRRm)
# 0.4239479
sd(l_vraiRRm)
# 0.007845985

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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, I1_prev = 0.3,
                 co_inf_para1 = 1, co_inf_para2 = -1)
  
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
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                     |
# |:----|:-------|:------------------------------------------|:------------------------------------------|
# |1000 |RegLog  |MCSE_bias          , 0.00264657206939816   |Bias             , 0.099309108007026       |
# |-    |-       |MCSE_var            , 0.000369571367405791 |Var                , 0.00700434371851846   |
# |-    |-       |MCSE_mse            , 0.000760064590618005 |Mse               , 0.0168596383079511     |
# |-    |-       |%Cov , 0.789                               |Précision_var       , -0.00340579269112104 |
  
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
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                     |
# |:----|:-------|:------------------------------------------|:------------------------------------------|
# |1000 |IPW     |MCSE_bias          , 0.00270552263609028   |Bias              , 0.0825316965966604     |
# |-    |-       |MCSE_var            , 0.000386930763215094 |Var                , 0.00731985273439689   |
# |-    |-       |MCSE_mse            , 0.000760064590618005 |Mse               , 0.0168596383079511     |
# |-    |-       |%Cov, 0.85                                 |Précision_var       , -0.00281724227395977 |

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
# |n    |Methode  |Erreur de Monte Carlo                      |Autres                                     |
# |:----|:--------|:------------------------------------------|:------------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00289581023606088   |Bias              , 0.0678870838330913     |
# |-    |-        |MCSE_var            , 0.000490830785858227 |Var                , 0.00838571692327496   |
# |-    |-        |MCSE_mse            , 0.000707781915721007 |Mse               , 0.0129859873577128     |
# |-    |-        |%Cov , 0.947                               |Précision_var       , -0.00381030944386027 |

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
# |1000 |TNDDR_Lasso |MCSE_bias          , 0.00509663874242636 |Bias             , 0.108386248161959   |
# |-    |-           |MCSE_var           , 0.00178930637904432 |Var               , 0.0259757264708014 |
# |-    |-           |MCSE_mse           , 0.00181213306696707 |Mse               , 0.0376973295349564 |
# |-    |-           |%Cov , 0.506                             |Précision_var    , 0.259092060773553   |

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
# |1000 |TNDDR_Mars |MCSE_bias         , 0.0052798559563943   |Bias            , 0.04118294201346     |
# |-    |-          |MCSE_var           , 0.00182424624794171 |Var               , 0.0278768789202723 |
# |-    |-          |MCSE_mse          , 0.0018045175714018   |Mse               , 0.0295450367542361 |
# |-    |-          |%Cov , 0.961                             |Précision_var  , 0.8611806884107       |

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
# |1000 |TNDDR_RN |MCSE_bias          , 0.00905372840210406 |Bias               , -0.0569631024018428 |
# |-    |-        |MCSE_var           , 0.00290789018373685 |Var               , 0.0819699979790657   |
# |-    |-        |MCSE_mse           , 0.00291069335129148 |Mse               , 0.0851328230163295   |
# |-    |-        |%Cov , 0.912                             |Précision_var    , 0.496274150298166     |

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
# |1000 |TNDDR_GLM |MCSE_bias         , 0.0122455467978932   |Bias              , -0.200460867029602 |
# |-    |-         |MCSE_var           , 0.00663990119311776 |Var              , 0.149953416379391   |
# |-    |-         |MCSE_mse           , 0.00242490558831578 |Mse              , 0.189988022173272   |
# |-    |-         |%Cov, 1                                  |Précision_var    , -11.1571025298495   |

################################################################################
################################################################################
# Pour I1_prev = 50% : co_inf_para1 = 1.5, co_inf_para2 = -1.5
# Proportion de co-infection ~30%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], I1_prev = 0.5,
                      co_inf_para1 = 1.5, co_inf_para2 = -1.5)
  
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
# 0.4214469
sd(l_vraiRRc)
# 0.008242331

l_vraiRRm

mean(l_vraiRRm)
# 0.4307477
sd(l_vraiRRm)
# 0.008206604

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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, I1_prev = 0.5,
                 co_inf_para1 = 1.5, co_inf_para2 = -1.5)
  
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
# |1000 |RegLog  |MCSE_bias          , 0.00273221254378159   |Bias             , 0.128218515635666     |
# |-    |-       |MCSE_var            , 0.000338596404089793 |Var                , 0.00746498538439746 |
# |-    |-       |MCSE_mse            , 0.000855421079617742 |Mse               , 0.0238975081508266   |
# |-    |-       |%Cov , 0.596                               |Précision_var      , 0.00813363371670661 |

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
# |1000 |IPW     |MCSE_bias          , 0.00272317326312777   |Bias             , 0.109212855829551     |
# |-    |-       |MCSE_var           , 0.00034314151424524   |Var                , 0.00741567262101392 |
# |-    |-       |MCSE_mse            , 0.000855421079617742 |Mse               , 0.0238975081508266   |
# |-    |-       |%Cov , 0.705                               |Précision_var      , 0.00738823760476939 |

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
# |1000 |TNDDR_RF |MCSE_bias          , 0.00286075909721386   |Bias              , 0.0938633917358341   |
# |-    |-        |MCSE_var            , 0.000469183967868548 |Var                , 0.00818394261229187 |
# |-    |-        |MCSE_mse            , 0.000766149432576603 |Mse               , 0.0169860949778342   |
# |-    |-        |%Cov , 0.882                               |Précision_var     , 0.0912015905663806   |

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
# |1000 |TNDDR_Lasso |MCSE_bias         , 0.0029602413718813     |Bias             , 0.125245272856772     |
# |-    |-           |MCSE_var            , 0.000528958947264412 |Var                , 0.00876302897979769 |
# |-    |-           |MCSE_mse            , 0.000915699271588437 |Mse               , 0.0244406443237852   |
# |-    |-           |%Cov , 0.369                               |Précision_var    , 0.158963270139735     |

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
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00566199333251944 |Bias              , 0.0652315618986083 |
# |-    |-          |MCSE_var           , 0.00209276156484285 |Var               , 0.0320581684974946 |
# |-    |-          |MCSE_mse           , 0.00202696169741705 |Mse              , 0.036281266996729   |
# |-    |-          |%Cov , 0.932                             |Précision_var   , 0.64696419272813     |

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
# |1000 |TNDDR_RN |MCSE_bias          , 0.00784367288285133 |Bias              , 0.0457540019633394 |
# |-    |-        |MCSE_var           , 0.00282948321391781 |Var               , 0.0615232042931773 |
# |-    |-        |MCSE_mse           , 0.00280448853248097 |Mse               , 0.0635551097845454 |
# |-    |-        |%Cov , 0.909                             |Précision_var   , 1.51216120352509     |

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
# |1000 |TNDDR_GLM |MCSE_bias         , 0.0126159441364819   |Bias              , -0.197716794440465 |
# |-    |-         |MCSE_var           , 0.00665875486061737 |Var              , 0.159162046454831   |
# |-    |-         |MCSE_mse           , 0.00229602993902565 |Mse             , 0.19809481521219     |
# |-    |-         |%Cov, 1                                  |Précision_var    , -9.39368551619764   |