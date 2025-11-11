# Chargement des librairies nécessaires

install.packages("geepack")
library(geepack)


# Initialiser des objets pour contenir les resultats
#    - Coefficient 
#    - Erreur-type du coefficient


resultats2 <- data.frame(matrix(ncol = 4, 
                               nrow = nsim))

colnames(resultats2) <- c("iteration",
                         "germe",
                         "RRm_IPW",
                         "VE") # Risque relatif marginal

# Faire une boucle de 1 a nombre de replications
# Generation des donnees TND (n = 1000)

seeds_list <- sample(1:1000000, size = nsim)

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para = 0)
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
  
  # Calcul des poids
  
  mod.denom <- glm(V ~ C,
                   
                   family = binomial(link = "logit"),  
                   data = TNDdat,
                   subset = (dat$Y == 0)) # Chez les témoins
  
  g1 <- predict(mod.denom, type = "response")

  # Analyse avec IPW
  # Enregistrer les resultats
  
    ## TND (infection symptomatique)
  
  RRm <- mean(TNDdat$Y*TNDdat$V/g1)/mean(TNDdat$Y*(1-TNDdat$V)/(1-g1))
  
  ## Numéro de l'itération et la valeur du germe
  resultats2[i, 1] <- i
  resultats2[i, 2] <- seeds_list[i]
  
  ## Risque relatif marginal
  resultats2[i, 3] <- RRm
  
  ## Efficacité vaccinale
  resultats2[i, 4] <- 1 - RRm
  
  if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  # DT : Ajout d'une ligne pour suivre l'avancement
  
}

# Analyser les resultats

#    - statistiques descriptives

summary(resultats2$RRm)
mean(resultats2$RRm)
sd(resultats2$RRm)

#    - biais, variance, moyenne de l'erreur-type, couverture des IC

## Chargements des librairies nécessaires

library(simhelpers)
library(dplyr)
library(tibble)
library(knitr)
library(dplyr)
library(kableExtra)

# Ajout de la colonne contenant la vraie valeur du paramètre

resultats2$vrai_param <- rep(0.3664181, nsim)

## Table 01 Résultats de l'étude de simulation (Scénario 01 :
##                                                    P_co-infections = 0.00001)

Tab01 <- data.frame(matrix(ncol = 3, 
                           nrow = 4))

colnames(Tab01) <- c("n",
                     "Methode", 
                     "IPW")

Tab01$n <- c("1000", "-", "-", "-")
Tab01$Methode <- c("MCSE_bias", "MCSE_var", "MCSE_mse", "%Cov")

estimations <- resultats2$RRm
vraie_valeur <- 0.3664181 # moyenne(autres_vraies_valeurs)

# Utiliser la fonction "calc_absolute" pour calculer les différentes mesures de performance

help("calc_absolute")

#Calculates absolute bias, variance, mean squared error (mse) and root mean squared error (rmse).
#The function also calculates the associated Monte Carlo standard errors.

T <- calc_absolute(resultats2, RRm, vrai_param, criteria = c("bias", "stddev", "rmse"))
kable(T, digits = 3)

### Calcul du biais

MCSE_biais <- calc_absolute(resultats2, RRm, vrai_param, criteria = "bias")
Tab01$IPW[1] <- MCSE_biais[3]

MCSE_biais <- sqrt(sum((estimations - mean(estimations))^2) / (nsim*(nsim - 1)))
Tab01$IPW[1] <- MCSE_biais

### Calcul de la variance

MCSE_var <- calc_absolute(resultats2, RRm, vrai_param, criteria = "var")
Tab01$IPW[2] <- MCSE_var[3]

### Calcul de "Monte Carlo Mean Squared Error" (MC MSE)

MCSE_MSE <- calc_absolute(resultats2, RRm, vrai_param, criteria = "mse")
Tab01$IPW[3] <- MCSE_MSE[3]

################################################################################
############################ m-estimateurs #####################################

dat <- datagen(seed = 123456, ssize = 1000, co_inf_para = 0)
TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)

library(geex)

# Estimation de la variance

#### Estimation du RRm avec GEEx ####

geex_ef <- function(data){
  
  Y <- data$Y
  V <- data$V
  C <- data$C
  
  function(theta){
    
    alpha <- theta[1:2]  
    
    pscore <- (Y == 0)*plogis(alpha[1] + alpha[2]*C)
    # pscore seulement chez les témoins
    
    # Equations d'estimation
    
    eq_1 <- (Y == 0)*((V - 1)*pscore - V)                  
    eq_2 <- (Y == 0)*((V - 1)*pscore - V)*C

    eq_3 <- (Y*V/plogis(alpha[1] + alpha[2]*C)) - theta[3]
    eq_4 <- (Y*(1 - V)/(1 - plogis(alpha[1] + alpha[2]*C))) - theta[4]
    
    return(c(eq_1, eq_2, eq_3, eq_4))
  }
}

mestr <- m_estimate(estFUN = geex_ef,                                       
                    data = TNDdat,                                             
                    root_control = setup_root_control(start = c(0, 0, 0.5, 0.5))) # Mêmes valeurs que dans le code que vous m'avez envoyé

beta_geex <- roots(mestr)             
se_geex <- sqrt(diag(vcov(mestr))) # vcov(mestr) : Matrice de variance-covariance

## Comparaison

beta_geex[5:6]

# Intervalles de confiance

bor_Inf <- beta_geex - 1.96 * se_geex
bor_Sup <- beta_geex + 1.96 * se_geex


IC <- cbind( # Intervalle de confiance
  
  Estimate = beta_geex,
  Std.Error = se_geex,
  Lower_CI = bor_Inf,
  Upper_CI = bor_Sup
)

print(confidence_intervals)
