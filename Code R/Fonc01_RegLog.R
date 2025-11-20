
# Chargements des librairies nécessaires

library(dplyr)

# Initialiser des objets pour contenir les resultats

resultats <- data.frame(matrix(ncol = 6, 
                               nrow = nsim))

colnames(resultats) <- c("coe_reg", 
                         "err_reg",
                         "RRc", # Risque relatif conditionnel
                         "est_VE")

##    - biais, variance, moyenne de l'erreur-type, couverture des IC

resultats$vrai_param <- rep(0.6507249, nsim) # vraie valeur du paramètre

Tab01 <- data.frame(matrix(ncol = 3, 
                           nrow = 4))

colnames(Tab01) <- c("n",
                     "Methode", 
                     "Regression_logistique")

Tab01$Methode <- c("MCSE_bias", "MCSE_var", "MCSE_mse", "%Cov")
# Définir la fonction pour la régression logistique

RegLog <- function(dat){
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
 
  # Analyse avec regression logistique
  
  fit.TND <- glm(Y ~ V + C,
                 
                 family = binomial(link = "logit"),  
                 # DT: Pour les donnees TND, on utilise logit
                 data = TNDdat) 
  
  # Enregistrer les resultats
  
  resultats.TND <- summary(fit.TND)
  
  ## Coefficient de la regression logistique
  ## Erreur-type du coefficient
  ## RRc
  ## ^VE
  
  l <- list(resultats.TND$coefficients[2], resultats.TND$coefficients[5],
            exp(resultats.TND$coefficients[2]), 1 - exp(resultats.TND$coefficients[2]))  
  

  return(l)
  
}

