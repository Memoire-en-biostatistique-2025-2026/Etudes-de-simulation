# 1. Calculer les vraies valeurs 
# 2. Parametres de la simulation
#    - Nombre de replications 
#    - Liste des germes
# 3. Initialiser des objets pour contenir les resultats
#    - Coefficient de la regression logistique
#    - Erreur-type du coefficient
# 4. Faire une boucle de 1 a nombre de replications
# 5. Generation des donnees TND 
# 6. Analyse avec regression logistique
# 7. Enregistrer les resultats
# 8. Analyser les resultats
#    - statistiques descriptives
#    - biais, variance, moyenne de l'erreur-type, couverture des IC
################################################################################
################################################################################

# 2. Parametres de la simulation
#    - Nombre de replications (10)
#    - Liste de 10 germes

## Regression logistique ##

## Création d'une liste de germes pour assurer la reproductibilité des résultats



set.seed(1) # Pour avoir toujours les mêmes germes

nsim <- 10 # Fixer le nombre de réplications

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nrep)
vrai.EV.autres <- rep(NA, nrep)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], popsize = 10000000)
  
  dat0 <- data.frame(C = dat$C, V = 0, Y = dat$I2_0*dat$W2_0*dat$H_0)
  
  dat1 <- data.frame(C = dat$C, V = 1, Y = dat$I2_1*dat$W2_1*dat$H_1) 
  
  dat_complet <- rbind(dat0, dat1)
  
  vraiRRc <- glm(Y ~ V + C, family = binomial(link = "log"), data = dat_complet)
  
  l_vraiRRc[i] <- 1 - exp(coef(vraiRRc)[2])
  
  vraiRRm <- mean(dat1$Y)/mean(dat0$Y)
  
  vrai.EV.autres[i] <- 1 - vraiRRm
  
  print(data.frame(Sys.time(), i))
  
}

l_vraiRRc
# [1] 0.3597205 0.3618331 0.3631053 0.3631609 0.3640605 0.3611498 0.3627725 0.3612469 0.3623667
# [10] 0.3632017
mean(l_vraiRRc)
# [1] 0.3622618
sd(l_vraiRRc)
# [1] 0.001283177

# 3. Initialiser des objets pour contenir les resultats
#    - Coefficient de la regression logistique
#    - Erreur-type du coefficient


resultats <- data.frame(matrix(ncol = 5, 
                            nrow = nsim))

colnames(resultats) <- c("iteration",
                         "germe",
                         "RRc", 
                         "err_reg",
                         "VE")

# 4. Faire une boucle de 1 a nombre de replications
# 5. Generation des donnees TND (n = 1000)

nsim <- 100 # Fixer le nombre de réplications

seeds_list <- sample(1:1000000, size = nsim)

for (i in 1:nsim) {
  
    dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para = 0)
    TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
    
# 6. Analyse avec regression logistique
# 7. Enregistrer les resultats
    
    ## TND (infection symptomatique)
    
    fit.TND <- glm(Y ~ V + C,
                       
                       family = binomial(link = "logit"),  
                        # DT: Pour les donnees TND, on utilise logit
                       data = TNDdat) 
    
    resultats.TND <- summary(fit.TND)
    
    ## Numéro de l'itération et la valeur du germe
    resultats[i, 1] <- i
    resultats[i, 2] <- seeds_list[i]
    
    ## Coefficient de la regression logistique
    resultats[i, 3] <- exp(resultats.TND$coefficients[2])

    ## Erreur-type du coefficient
    resultats[i, 4] <- resultats.TND$coefficients[5]
    
    resultats[i, 5] <- 1 - exp(resultats.TND$coefficients[2])
    
    if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
    
      # DT : Ajout d'une ligne pour suivre l'avancement

}

# 8. Analyser les resultats

#    - statistiques descriptives

summary(resultats$RRc)
mean(resultats$RRc)
sd(resultats$RRc)

#    - biais, variance, moyenne de l'erreur-type, couverture des IC

## Chargements des librairies nécessaires

library(simhelpers)
library(dplyr)
library(tibble)
library(knitr)
library(dplyr)
library(kableExtra)

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats$vrai_param <- rep(0.3622618, nsim)

## Table 01 Résultats de l'étude de simulation (Scénario 01 :
##                                                    P_co-infections = 0.00001)

Tab01 <- data.frame(matrix(ncol = 3, 
                               nrow = 4))

colnames(Tab01) <- c("n",
                      "Methode", 
                      "Regression_logistique")

Tab01$n <- c("1000", "-", "-", "-")
Tab01$Methode <- c("MCSE_bias", "MCSE_var", "MCSE_mse", "%Cov")

estimations <- resultats$RRc
vraie_valeur <- 0.3622618 # moyenne(vraies_valeurs)

# Utiliser la fonction "calc_absolute" pour calculer les différentes mesures de performance

help("calc_absolute")

#Calculates absolute bias, variance, mean squared error (mse) and root mean squared error (rmse).
#The function also calculates the associated Monte Carlo standard errors.

T <- calc_absolute(resultats, RRc, vrai_param, criteria = c("bias", "stddev", "rmse"))
kable(T, digits = 3)

### Calcul du biais

MCSE_biais <- calc_absolute(resultats, RRc, vrai_param, criteria = "bias")
Tab01$Regression_logistique[1] <- MCSE_biais[3]

MCSE_biais <- sqrt(sum((estimations - mean(estimations))^2) / (nsim*(nsim - 1)))
Tab01$Regression_logistique[1] <- MCSE_biais

### Calcul de la variance

MCSE_var <- calc_absolute(resultats, RRc, vrai_param, criteria = "var")
Tab01$Regression_logistique[2] <- MCSE_var[3]

### Calcul de "Monte Carlo Mean Squared Error" (MC MSE)

MCSE_MSE <- calc_absolute(resultats, RRc, vrai_param, criteria = "mse")
Tab01$Regression_logistique[3] <- MCSE_MSE[3]

### Calcul de la "couverture" (coverage)

help("calc_coverage")

# Calcul des bornes de l'intervalle de confiance

resultats$lim_inf <- 1 - exp(resultats[, 3] + 1.96*resultats[, 4])
resultats$lim_sup <- 1 - exp(resultats[, 3] - 1.96*resultats[, 4])

mean(resultats$lim_inf < resultats$vrai_param & resultats$lim_sup > resultats$vrai_param)
  
coverage <- calc_coverage(resultats, lim_inf, lim_sup, vrai_param)
Tab01$Regression_logistique[4] <- coverage[2]

kable(Tab01)
