
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

# Calculate true values of cRR and mRR 

set.seed(1) # To ensure reproducibility

nsim <- 10 # Number of replicas

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i],
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
# 0.433797
sd(l_vraiRRc)
# 0.01140663

l_vraiRRm

mean(l_vraiRRm)
# 0.442533
sd(l_vraiRRm)
# 0.01140144

################################################################################
# Initialize objects to store results

Tab01 <- data.frame(n = c("5000", "-", "-", "-"))

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

# Initialize objects to store results for replications with very different results

replicats <- data.frame(matrix(ncol = 9, 
                               nrow = 3*nsim))

colnames(replicats) <- c("n°Replication", # Numéro de la réplicatoon
                         "Méthode", # Méthode d'estimation
                         "mu1", "mu0",
                         "m0",
                         "g1", "g0",
                         "w1", "w0") 


methode <- list(RandomForest, Lasso, Mars, RN, PM)

r <- 1

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 5000, co_inf_para1 = 8, co_inf_para2 = -8)
  
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
  
  if((resultats3$RRm_Lasso[i] - resultats3$RRm_Mars[i]) > 0.1){
    
    replicats[r,c(1, 2)] <- list(i, "Lasso")
    replicats[r, -c(1, 2)] <- Lasso(dat)
    r <- r + 1 
    replicats[r,c(1, 2)] <- list(i, "GLM")
    replicats[r, -c(1, 2)] <- PM(dat)
    r <- r + 1
    replicats[r,c(1, 2)] <- list(i, "Mars")
    replicats[r, -c(1, 2)] <- Mars(dat)
    r <- r + 1
    
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
# |n    |Methode |Erreur de Monte Carlo                      |Autres                                       |
# |:----|:-------|:------------------------------------------|:--------------------------------------------|
# |5000 |RegLog  |MCSE_bias          , 0.00108431918588206   |Bias                 , -0.000851296641909416 |
# |-    |-       |MCSE_var            , 5.45949418635126e-05 |Var                , 0.00117574809687194     |
# |-    |-       |MCSE_mse            , 5.44101253214418e-05 |Mse               , 0.0011752970547476       |
# |-    |-       |%Cov , 0.923                               |Précision_var     , 0.0062566663507647       |

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
# |5000 |IPW     |MCSE_bias          , 0.00115433684714761   |Bias               , -0.0025517212516738 |
# |-    |-       |MCSE_var            , 6.25208607121234e-05 |Var                , 0.00133249355668268 |
# |-    |-       |MCSE_mse            , 5.44101253214418e-05 |Mse               , 0.0011752970547476   |
# |-    |-       |%Cov , 0.911                               |Précision_var      , 0.00760524634280019 |    |-       |%Cov , 0.911                               |Précision_var      , 0.00760524634280019 |

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
# |5000 |TNDDR_RF |MCSE_bias          , 0.00117487695532648   |Bias                , -0.00675788529298121 |
# |-    |-        |MCSE_var            , 6.78102694165151e-05 |Var                , 0.00138033586015723   |
# |-    |-        |MCSE_mse            , 6.63329863403153e-05 |Mse                , 0.00142462453793016   |
# |-    |-        |%Cov , 0.891                               |Précision_var      , 0.00993714826852501   |

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
# |5000 |TNDDR_Lasso |MCSE_bias          , 0.00117212478677664   |Bias               , -0.0029331079850875 |
# |-    |-           |MCSE_var            , 6.68833767689461e-05 |Var                , 0.00137387651577618 |
# |-    |-           |MCSE_mse            , 6.59562185055415e-05 |Mse                , 0.00138110576171259 |
# |-    |-           |%Cov, 0.81                                 |Précision_var     , 0.0267239136581471   |

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
# |n    |Methode    |Erreur de Monte Carlo                      |Autres                                   |
# |:----|:----------|:------------------------------------------|:----------------------------------------|
# |5000 |TNDDR_Mars |MCSE_bias          , 0.00240036895814377   |Bias               , -0.0130665749811152 |
# |-    |-          |MCSE_var            , 0.000838485381456876 |Var                , 0.00576177113522019 |
# |-    |-          |MCSE_mse            , 0.000857953551471391 |Mse                , 0.00592674474582207 |
# |-    |-          |%Cov , 0.917                               |Précision_var    , 0.545987782999199     |

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
# |n    |Methode  |Erreur de Monte Carlo                      |Autres                                     |
# |:----|:--------|:------------------------------------------|:------------------------------------------|
# |5000 |TNDDR_RN |MCSE_bias          , 0.00136982059854548   |Bias                , -0.00706337957076952 |
# |-    |-        |MCSE_var            , 0.000243174790190174 |Var                , 0.00187640847219949   |
# |-    |-        |MCSE_mse            , 0.000251100049000652 |Mse                , 0.00192442339468806   |
# |-    |-        |%Cov , 0.904                               |Précision_var   , 0.13435367339021         |

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
# |5000 |TNDDR_GLM |MCSE_bias          , 0.00116479868635591   |Bias                , -0.00347919893509951 |
# |-    |-         |MCSE_var            , 6.50844275131659e-05 |Var                , 0.00135675597973645   |
# |-    |-         |MCSE_mse           , 6.4216453404653e-05   |Mse                , 0.00136750404898671   |
# |-    |-         |%Cov, 0.91                                 |Précision_var      , 0.00649362492378254   |