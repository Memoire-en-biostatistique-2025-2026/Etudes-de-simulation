
# Chargements des librairies nécessaires

library(dplyr)


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

