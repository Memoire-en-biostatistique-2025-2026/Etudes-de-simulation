
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
# Scénario 05 (co-infections ~ 40% in the TND sample): Sous-scénario 01 : I2_prev = 10%, 15% et 30% 
#                                                                         I1_prev = 15% et 
#                                                                         CV = 33%
# For I2_prev = 10% : co_inf_para1 = 3, co_inf_para2 = 5
# For I2_prev = 15% : co_inf_para1 = 3, co_inf_para2 = 5
# For I2_prev = 30% : co_inf_para1 = 1, co_inf_para2 = 0.5

################################################################################

################################################################################
################################################################################
# For I2_prev = 10% : co_inf_para1 = 3, co_inf_para2 = 5
# Proportion of co-infections ~ 40%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], I2_prev = 0.1,
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
# 0.4128448
sd(l_vraiRRc)
# 0.008443783

l_vraiRRm

mean(l_vraiRRm)
# 0.4149135
sd(l_vraiRRm)
# 0.008455331

################################################################################
# Initialize objects to store results

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

################################################################################
########################## Analyse des résultats ###############################

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
                 co_inf_para1 = 3, co_inf_para2 = 5)
  
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

############################# Logistic regression ##############################

Tab01$Methode <- c("RegLog", "-", "-","-")

##    - Descriptive statistics

summary(resultats$RRc)
mean(resultats$RRc)
sd(resultats$RRc)

sd(resultats$coe_reg)
mean(resultats$err_reg)

help("calc_absolute") # calculate various performance metrics

# Add a column containing the true value of the parameter

resultats$vrai_param <- rep(mean(l_vraiRRc), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats, RRc, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats, RRc, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats, RRc, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calculating the confidence interval bounds

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
  list(name = "Précision_var", value = sd(resultats$coe_reg) - mean(resultats$err_reg)) # See how accurately the variance is estimated
  
)

kable(Tab01)
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-------|:------------------------------------------|:----------------------------------------|
# |1000 |RegLog  |MCSE_bias          , 0.00259780044976044   |Bias             , 0.116253262677273     |
# |-    |-       |MCSE_var            , 0.000303026770348753 |Var                , 0.00674856717677554 |
# |-    |-       |MCSE_mse            , 0.000747373296352107 |Mse               , 0.0202566396927097   |
# |-    |-       |%Cov , 0.611                               |Précision_var     , 0.0093910277046671   |

################################ IPW ###########################################

Tab01$Methode <- c("IPW", "-", "-","-")

##    - Descriptive statistics

summary(resultats2$RRm)
mean(resultats2$RRm)
sd(resultats2$RRm)

mean(sqrt(resultats2$var_log_RRm))
sd(log(resultats2$RRm))

help("calc_absolute") # calculate various performance metrics

# Add a column containing the true value of the parameter

resultats2$vrai_param <- rep(mean(l_vraiRRm), nsim) 

### MCSE_biais

MCSE_biais <- calc_absolute(resultats2, RRm, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats2, RRm, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_MSE <- calc_absolute(resultats2, RRm, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

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
  list(name = "Précision_var", value = sd(log(resultats2$RRm)) - mean(sqrt(resultats2$var_log_RRm)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-------|:------------------------------------------|:----------------------------------------|
# |1000 |IPW     |MCSE_bias          , 0.00255869289826532   |Bias             , 0.108370842857661     |
# |-    |-       |MCSE_var            , 0.000298365775160736 |Var                , 0.00654690934763339 |
# |-    |-       |MCSE_mse            , 0.000747373296352107 |Mse               , 0.0202566396927097   |
# |-    |-       |%Cov , 0.656                               |Précision_var     , 0.0078199990499061   |

################################## TNDDR #######################################

Tab01$Methode <- c("TNDDR_RF", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_RF)
mean(resultats3$RRm_RF)
sd(resultats3$RRm_RF)

mean(sqrt(resultats3$`var_log_RRm-RF`))
sd(log(resultats3$RRm_RF))

help("calc_absolute") # calculate various performance metrics

# Add a column containing the true value of the parameter

resultats3$vrai_param <- rep(mean(l_vraiRRm), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RF)) - mean(sqrt(resultats3$`var_log_RRm-RF`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:--------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00271083653990644   |Bias              , 0.0965544041220106   |
# |-    |-        |MCSE_var            , 0.000340719803817594 |Var                , 0.00734863474609192 |
# |-    |-        |MCSE_mse            , 0.000710567510575639 |Mse               , 0.0166640390667024   |
# |-    |-        |%Cov , 0.826                               |Précision_var      , 0.00707649472964758 |

################################################################################

# Lasso regression

Tab01$Methode <- c("TNDDR_Lasso", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_Lasso)
mean(resultats3$RRm_Lasso)
sd(resultats3$RRm_Lasso)

mean(sqrt(resultats3$`var_log_RRm-Lasso`))
sd(log(resultats3$RRm_Lasso))

help("calc_absolute") # calculate various performance metrics

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Lasso)) - mean(sqrt(resultats3$`var_log_RRm-Lasso`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode     |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-----------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias          , 0.00258798767194816   |Bias             , 0.104247341388901     |
# |-    |-           |MCSE_var           , 0.00031255765367242   |Var                , 0.00669768019015567 |
# |-    |-           |MCSE_mse            , 0.000696390809911209 |Mse               , 0.0175584906966196   |
# |-    |-           |%Cov , 0.679                               |Précision_var     , 0.0285688137569019   |

################################################################################

# MARS

Tab01$Methode <- c("TNDDR_Mars", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_Mars)
mean(resultats3$RRm_Mars)
sd(resultats3$RRm_Mars)

mean(sqrt(resultats3$`var_log_RRm-Mars`))
sd(log(resultats3$RRm_Mars))

help("calc_absolute") # calculate various performance metrics

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Mars)) - mean(sqrt(resultats3$`var_log_RRm-Mars`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode    |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:----------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00490205486872357 |Bias              , 0.0854137536787844 |
# |-    |-          |MCSE_var           , 0.00177655752205221 |Var               , 0.0240301419359765 |
# |-    |-          |MCSE_mse           , 0.00181821428506868 |Mse               , 0.0313016211115405 |
# |-    |-          |%Cov , 0.889                             |Précision_var    , 0.400958391012882   |

################################################################################

# Neural networks

Tab01$Methode <- c("TNDDR_RN", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_RN)
mean(resultats3$RRm_RN)
sd(resultats3$RRm_RN)

mean(sqrt(resultats3$`var_log_RRm-RN`))
sd(log(resultats3$RRm_RN))

help("calc_absolute") # calculate various performance metrics

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RN)) - mean(sqrt(resultats3$`var_log_RRm-RN`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:--------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_RN |MCSE_bias          , 0.00591506943749588 |Bias              , 0.0656971800293445 |
# |-    |-        |MCSE_var          , 0.0021765837184794   |Var               , 0.0349880464503978 |
# |-    |-        |MCSE_mse           , 0.00213739394918378 |Mse               , 0.0392691778677555 |
# |-    |-        |%Cov, 0.91                               |Précision_var    , 0.964595424470064   |

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
# |1000 |TNDDR_GLM |MCSE_bias          , 0.00259051107160595   |Bias             , 0.105906711636692     |
# |-    |-         |MCSE_var            , 0.000310091143127994 |Var                , 0.00671074761211299 |
# |-    |-         |MCSE_mse            , 0.000699174362314067 |Mse               , 0.0179202684341983   |
# |-    |-         |%Cov , 0.797                               |Précision_var      , 0.00592233956796778 |

################################################################################
################################################################################
# For I2_prev = 15% : co_inf_para1 = 3, co_inf_para2 = 5
# Proportion of co-infections ~ 30%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], I2_prev = 0.15,
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
# 0.416099
sd(l_vraiRRc)
# 0.01066708

l_vraiRRm

mean(l_vraiRRm)
# 0.4191034
sd(l_vraiRRm)
# 0.01068322

################################################################################
# Initialize objects to store results

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

################################################################################
########################## Analyse des résultats ###############################

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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, I2_prev = 0.15,
                 co_inf_para1 = 3, co_inf_para2 = 5)
  
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

############################# Logistic regression ##############################

Tab01$Methode <- c("RegLog", "-", "-","-")

##    - Descriptive statistics

summary(resultats$RRc)
mean(resultats$RRc)
sd(resultats$RRc)

sd(resultats$coe_reg)
mean(resultats$err_reg)

help("calc_absolute") # calculate various performance metrics

# Add a column containing the true value of the parameter

resultats$vrai_param <- rep(mean(l_vraiRRc), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats, RRc, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats, RRc, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats, RRc, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calculating the confidence interval bounds

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
  list(name = "Précision_var", value = sd(resultats$coe_reg) - mean(resultats$err_reg)) # See how accurately the variance is estimated
  
)

kable(Tab01)
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-------|:------------------------------------------|:----------------------------------------|
# |1000 |RegLog  |MCSE_bias          , 0.00262903452818251   |Bias             , 0.119426707648531     |
# |-    |-       |MCSE_var            , 0.000303598013532188 |Var                , 0.00691182255037583 |
# |-    |-       |MCSE_mse            , 0.000753079917884466 |Mse               , 0.0211676492275932   |
# |-    |-       |%Cov, 0.61                                 |Précision_var      , 0.00852016698688271 |
  
################################ IPW ###########################################

Tab01$Methode <- c("IPW", "-", "-","-")

##    - Descriptive statistics

summary(resultats2$RRm)
mean(resultats2$RRm)
sd(resultats2$RRm)

mean(sqrt(resultats2$var_log_RRm))
sd(log(resultats2$RRm))

help("calc_absolute") # calculate various performance metrics

# Add a column containing the true value of the parameter

resultats2$vrai_param <- rep(mean(l_vraiRRm), nsim) 

### MCSE_biais

MCSE_biais <- calc_absolute(resultats2, RRm, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats2, RRm, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_MSE <- calc_absolute(resultats2, RRm, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

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
  list(name = "Précision_var", value = sd(log(resultats2$RRm)) - mean(sqrt(resultats2$var_log_RRm)) # See how accurately the variance is estimated
       
  )
  
)


kable(Tab01)

# |n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-------|:------------------------------------------|:----------------------------------------|
# |1000 |IPW     |MCSE_bias          , 0.00257338085530473   |Bias             , 0.107480575053143     |
# |-    |-       |MCSE_var            , 0.000288960453371646 |Var                , 0.00662228902644892 |
# |-    |-       |MCSE_mse            , 0.000753079917884466 |Mse               , 0.0211676492275932   |
# |-    |-       |%Cov , 0.671                               |Précision_var      , 0.00572956513463066 |

################################## TNDDR #######################################

Tab01$Methode <- c("TNDDR_RF", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_RF)
mean(resultats3$RRm_RF)
sd(resultats3$RRm_RF)

mean(sqrt(resultats3$`var_log_RRm-RF`))
sd(log(resultats3$RRm_RF))

help("calc_absolute") # calculate various performance metrics

# Add a column containing the true value of the parameter

resultats3$vrai_param <- rep(mean(l_vraiRRm), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RF)) - mean(sqrt(resultats3$`var_log_RRm-RF`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                      |Autres                                       |
# |:----|:--------|:------------------------------------------|:--------------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00264876510622258   |Bias              , 0.0917729771328847       |
# |-    |-        |MCSE_var            , 0.000318200076554763 |Var                , 0.00701595658794229     |
# |-    |-        |MCSE_mse            , 0.000629318844248589 |Mse               , 0.0154312199631873       |
# |-    |-        |%Cov , 0.867                               |Précision_var        , -0.000210419364979869 |

################################################################################

# Lasso regression

Tab01$Methode <- c("TNDDR_Lasso", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_Lasso)
mean(resultats3$RRm_Lasso)
sd(resultats3$RRm_Lasso)

mean(sqrt(resultats3$`var_log_RRm-Lasso`))
sd(log(resultats3$RRm_Lasso))

help("calc_absolute") # calculate various performance metrics

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Lasso)) - mean(sqrt(resultats3$`var_log_RRm-Lasso`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode     |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-----------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias          , 0.00258974452740089   |Bias              , 0.0987122190682714   |
# |-    |-           |MCSE_var            , 0.000292313309984534 |Var                , 0.00670677671720287 |
# |-    |-           |MCSE_mse            , 0.000635385131268093 |Mse               , 0.0164441721338681   |
# |-    |-           |%Cov , 0.693                               |Précision_var    , 0.032292546985225     |

################################################################################

# MARS

Tab01$Methode <- c("TNDDR_Mars", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_Mars)
mean(resultats3$RRm_Mars)
sd(resultats3$RRm_Mars)

mean(sqrt(resultats3$`var_log_RRm-Mars`))
sd(log(resultats3$RRm_Mars))

help("calc_absolute") # calculate various performance metrics

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Mars)) - mean(sqrt(resultats3$`var_log_RRm-Mars`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode    |Erreur de Monte Carlo                    |Autres                                   |
# |:----|:----------|:----------------------------------------|:----------------------------------------|
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00518699193752641 |Bias              , 0.0737461067122665   |
# |-    |-          |MCSE_var           , 0.00184208124803755 |Var              , 0.026904885359964     |
# |-    |-          |MCSE_mse           , 0.00182801875271719 |Mse              , 0.032316468729821     |
# |-    |-          |%Cov , 0.923                             |Précision_var      , -0.0456928988706264 |

################################################################################

# Neural networks

Tab01$Methode <- c("TNDDR_RN", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_RN)
mean(resultats3$RRm_RN)
sd(resultats3$RRm_RN)

mean(sqrt(resultats3$`var_log_RRm-RN`))
sd(log(resultats3$RRm_RN))

help("calc_absolute") # calculate various performance metrics

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RN)) - mean(sqrt(resultats3$`var_log_RRm-RN`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:--------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_RN |MCSE_bias          , 0.00667147722728547 |Bias              , 0.0290988963277247 |
# |-    |-        |MCSE_var           , 0.00229376779250128 |Var               , 0.0445086083941886 |
# |-    |-        |MCSE_mse           , 0.00220722465261325 |Mse               , 0.0453108455532861 |
# |-    |-        |%Cov , 0.925                             |Précision_var   , 1.11567955772913     |

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
# |1000 |TNDDR_GLM |MCSE_bias          , 0.00263254443018267   |Bias             , 0.102354496484101     |
# |-    |-         |MCSE_var            , 0.000305864981115521 |Var                , 0.00693029017688579 |
# |-    |-         |MCSE_mse            , 0.000674524744368665 |Mse               , 0.0173998028372227   |
# |-    |-         |%Cov , 0.838                               |Précision_var      , 0.00289313043869857 |

################################################################################
################################################################################
# For I2_prev = 30% : co_inf_para1 = 1, co_inf_para2 = 0.5
# Proportion of co-infections ~ 30%

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i],  I2_prev = 0.3,
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
# 0.4124463
sd(l_vraiRRc)
# 0.009077658

l_vraiRRm

mean(l_vraiRRm)
# 0.4184352
sd(l_vraiRRm)
# 0.009015059

################################################################################
# Initialize objects to store results

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

################################################################################
########################## Analyse des résultats ###############################

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
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, I2_prev = 0.3,
                 co_inf_para1 = 1, co_inf_para2 = 0.5)
  
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

############################# Logistic regression ##############################

Tab01$Methode <- c("RegLog", "-", "-","-")

##    - Descriptive statistics

summary(resultats$RRc)
mean(resultats$RRc)
sd(resultats$RRc)

sd(resultats$coe_reg)
mean(resultats$err_reg)

help("calc_absolute") # calculate various performance metrics

# Add a column containing the true value of the parameter

resultats$vrai_param <- rep(mean(l_vraiRRc), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats, RRc, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats, RRc, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats, RRc, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

# Calculating the confidence interval bounds

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
  list(name = "Précision_var", value = sd(resultats$coe_reg) - mean(resultats$err_reg)) # See how accurately the variance is estimated
  
)

kable(Tab01)
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-------|:------------------------------------------|:----------------------------------------|
# |1000 |RegLog  |MCSE_bias          , 0.00271404574140266   |Bias              , 0.0861863901408258   |
# |-    |-       |MCSE_var            , 0.000318121924859839 |Var                , 0.00736604428642589 |
# |-    |-       |MCSE_mse            , 0.000631637150772615 |Mse               , 0.0147867720876461   |
# |-    |-       |%Cov, 0.81                                 |Précision_var      , 0.00495838898220402 |

################################ IPW ###########################################

Tab01$Methode <- c("IPW", "-", "-","-")

##    - Descriptive statistics

summary(resultats2$RRm)
mean(resultats2$RRm)
sd(resultats2$RRm)

mean(sqrt(resultats2$var_log_RRm))
sd(log(resultats2$RRm))

help("calc_absolute") # calculate various performance metrics

# Add a column containing the true value of the parameter

resultats2$vrai_param <- rep(mean(l_vraiRRm), nsim) 

### MCSE_biais

MCSE_biais <- calc_absolute(resultats2, RRm, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats2, RRm, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_MSE <- calc_absolute(resultats2, RRm, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

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
  list(name = "Précision_var", value = sd(log(resultats2$RRm)) - mean(sqrt(resultats2$var_log_RRm)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-------|:------------------------------------------|:----------------------------------------|
# |1000 |IPW     |MCSE_bias         , 0.0027214569187437     |Bias              , 0.0691235305956436   |
# |-    |-       |MCSE_var            , 0.000340366267780975 |Var                , 0.00740632776057794 |
# |-    |-       |MCSE_mse            , 0.000631637150772615 |Mse               , 0.0147867720876461   |
# |-    |-       |%Cov , 0.881                               |Précision_var     , 0.0025036781698787   |

################################## TNDDR #######################################

Tab01$Methode <- c("TNDDR_RF", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_RF)
mean(resultats3$RRm_RF)
sd(resultats3$RRm_RF)

mean(sqrt(resultats3$`var_log_RRm-RF`))
sd(log(resultats3$RRm_RF))

help("calc_absolute") # calculate various performance metrics

# Add a column containing the true value of the parameter

resultats3$vrai_param <- rep(mean(l_vraiRRm), nsim)

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_RF, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RF)) - mean(sqrt(resultats3$`var_log_RRm-RF`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                      |Autres                                     |
# |:----|:--------|:------------------------------------------|:------------------------------------------|
# |1000 |TNDDR_RF |MCSE_bias          , 0.00284293809862495   |Bias              , 0.0556382038566711     |
# |-    |-        |MCSE_var            , 0.000368588949550095 |Var                , 0.00808229703261325   |
# |-    |-        |MCSE_mse            , 0.000544698315154747 |Mse               , 0.0111698244639771     |
# |-    |-        |%Cov , 0.955                               |Précision_var       , -0.00426981535766313 |

################################################################################

# Lasso regression

Tab01$Methode <- c("TNDDR_Lasso", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_Lasso)
mean(resultats3$RRm_Lasso)
sd(resultats3$RRm_Lasso)

mean(sqrt(resultats3$`var_log_RRm-Lasso`))
sd(log(resultats3$RRm_Lasso))

help("calc_absolute") # calculate various performance metrics

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_Lasso, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Lasso)) - mean(sqrt(resultats3$`var_log_RRm-Lasso`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode     |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:-----------|:------------------------------------------|:----------------------------------------|
# |1000 |TNDDR_Lasso |MCSE_bias          , 0.00276268876703406   |Bias              , 0.0568300230207686   |
# |-    |-           |MCSE_var            , 0.000359618528891898 |Var                , 0.00763244922349616 |
# |-    |-           |MCSE_mse            , 0.000537365383073602 |Mse               , 0.0108544682908137   |
# |-    |-           |%Cov , 0.863                               |Précision_var     , 0.0419909179114444   |

################################################################################

# MARS

Tab01$Methode <- c("TNDDR_Mars", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_Mars)
mean(resultats3$RRm_Mars)
sd(resultats3$RRm_Mars)

mean(sqrt(resultats3$`var_log_RRm-Mars`))
sd(log(resultats3$RRm_Mars))

help("calc_absolute") # calculate various performance metrics

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_Mars, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_Mars)) - mean(sqrt(resultats3$`var_log_RRm-Mars`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode    |Erreur de Monte Carlo                    |Autres                                 |
# |:----|:----------|:----------------------------------------|:--------------------------------------|
# |1000 |TNDDR_Mars |MCSE_bias          , 0.00568641103582396 |Bias              , 0.0224434965882367 |
# |-    |-          |MCSE_var           , 0.00196876468135847 |Var               , 0.0323352704683405 |
# |-    |-          |MCSE_mse           , 0.00197366276937866 |Mse               , 0.0328066457369784 |
# |-    |-          |%Cov , 0.956                             |Précision_var    , 0.777796825601235   |

################################################################################

# Neural networks

Tab01$Methode <- c("TNDDR_RN", "-", "-","-")

##    - Descriptive statistics

summary(resultats3$RRm_RN)
mean(resultats3$RRm_RN)
sd(resultats3$RRm_RN)

mean(sqrt(resultats3$`var_log_RRm-RN`))
sd(log(resultats3$RRm_RN))

help("calc_absolute") # calculate various performance metrics

### MCSE_biais

MCSE_biais <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "bias")

### MCSE_var

MCSE_var <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "var")

### MCSE_MSE 

MCSE_mse <- calc_absolute(resultats3, RRm_RN, vrai_param, criteria = "mse")

### Coverage

help("calc_coverage")

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
  list(name = "Précision_var", value = sd(log(resultats3$RRm_RN)) - mean(sqrt(resultats3$`var_log_RRm-RN`)) # See how accurately the variance is estimated
       
  )
  
)

kable(Tab01)
# |n    |Methode  |Erreur de Monte Carlo                    |Autres                                   |
# |:----|:--------|:----------------------------------------|:----------------------------------------|
# |1000 |TNDDR_RN |MCSE_bias          , 0.00884817513233954 |Bias               , -0.0494875454410889 |
# |-    |-        |MCSE_var           , 0.00295768492059093 |Var               , 0.0782902031725518   |
# |-    |-        |MCSE_mse           , 0.00289212880553986 |Mse               , 0.0806609301231631   |
# |-    |-        |%Cov , 0.936                             |Précision_var   , 1.08816036078053       |

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
# |n    |Methode   |Erreur de Monte Carlo                      |Autres                                     |
# |:----|:---------|:------------------------------------------|:------------------------------------------|
# |1000 |TNDDR_GLM |MCSE_bias          , 0.00280457558975386   |Bias              , 0.0598054156496454     |
# |-    |-         |MCSE_var           , 0.00035546071549717   |Var                , 0.00786564423864319   |
# |-    |-         |MCSE_mse            , 0.000556499697765266 |Mse               , 0.0114344663354314     |
# |-    |-         |%Cov , 0.959                               |Précision_var       , -0.00490511062209914 |