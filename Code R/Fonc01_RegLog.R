
# Chargements des librairies nécessaires

library(simhelpers)
library(dplyr)
library(tibble)
library(knitr)
library(dplyr)
library(kableExtra)

# Initialiser des objets pour contenir les resultats

resultats <- data.frame(matrix(ncol = 6, 
                               nrow = nsim))

colnames(resultats) <- c("iteration",
                         "germe",
                         "coe_reg", 
                         "err_reg",
                         "RRc", # Risuqe relatif conditionnel
                         "est_VE")

seeds_list <- sample(1:1000000, size = nsim)

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000)
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
  
  # Analyse avec regression logistique
  
  fit.TND <- glm(Y ~ V + C,
                 
                 family = binomial(link = "logit"),  
                 # DT: Pour les donnees TND, on utilise logit
                 data = TNDdat) 
  
  # Enregistrer les resultats
  
  resultats.TND <- summary(fit.TND)
  
  ## Numéro de l'itération et la valeur du germe
  resultats[i, 1] <- i
  resultats[i, 2] <- seeds_list[i]
  
  ## Coefficient de la regression logistique
  resultats[i, 3] <- resultats.TND$coefficients[2]
  
  ## Erreur-type du coefficient
  resultats[i, 4] <- resultats.TND$coefficients[5]
  
  ## RRc
  resultats[i, 5] <- exp(resultats.TND$coefficients[2])
  
  ## ^VE
  resultats[i, 6] <- 1 - exp(resultats.TND$coefficients[2])
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
  # DT : Ajout d'une ligne pour suivre l'avancement
  
}

# Analyser les resultats

##    - statistiques descriptives

summary(resultats$RRc)
mean(resultats$RRc)
sd(resultats$RRc)

##    - biais, variance, moyenne de l'erreur-type, couverture des IC

resultats$vrai_param <- rep(0.6507249, nsim) # vraie valeur du paramètre

Tab01 <- data.frame(matrix(ncol = 3, 
                           nrow = 4))

colnames(Tab01) <- c("n",
                     "Methode", 
                     "Regression_logistique")

Tab01$n <- c("1000", "-", "-", "-")
Tab01$Methode <- c("MCSE_bias", "MCSE_var", "MCSE_mse", "%Cov")

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats, RRc, vrai_param, criteria = "bias")
Tab01$Regression_logistique[1] <- MCSE_biais[3]

### MCSE_var

MCSE_var <- calc_absolute(resultats, RRc, vrai_param, criteria = "var")
Tab01$Regression_logistique[2] <- MCSE_var[3]

### MCSE_MSE 

MCSE_MSE <- calc_absolute(resultats, RRc, vrai_param, criteria = "mse")
Tab01$Regression_logistique[3] <- MCSE_MSE[3]

### Coverage

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

resultats$lim_inf <- exp(resultats[, 3] - 1.96*resultats[, 4]) 
resultats$lim_sup <- exp(resultats[, 3] + 1.96*resultats[, 4])

mean(resultats$lim_inf < resultats$vrai_param & resultats$lim_sup > resultats$vrai_param)

coverage <- calc_coverage(resultats, lim_inf, lim_sup, vrai_param)
Tab01$Regression_logistique[4] <- coverage[2] # %Cov = 78%

kable(Tab01)


