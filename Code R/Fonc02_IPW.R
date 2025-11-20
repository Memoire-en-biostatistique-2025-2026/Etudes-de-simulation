
# Chargement des librairies nécessaires

library(dplyr)
library(geepack)
library(geex)

# Initialiser des objets pour contenir les resultats

resultats2 <- data.frame(matrix(ncol = 5, 
                                nrow = nsim))

colnames(resultats2) <- c("RRm_IPW",# Risque relatif marginal
                          "VE",
                          "var_RRm",# Variance du risque relatif marginal
                          "IC_inf", # Borne inférieure de l'intervalle de confiance
                          "IC_sup") # Sa borne supérieure


IPW <- function(dat){
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
  
  # Calcul des poids
  
  mod.denom <- glm(V ~ C,
                   
                   family = binomial(link = "logit"),  
                   data = TNDdat,
                   subset = (TNDdat$Y == 0)) # Chez les témoins
  
  g1 <- predict(mod.denom, type = "response")
  
  # Analyse avec IPW
  # Enregistrer les resultats
  
  ## TND (infection symptomatique)
  
  RRm <- mean(TNDdat$Y*TNDdat$V/g1)/mean(TNDdat$Y*(1-TNDdat$V)/(1-g1))
   
  # Risque relatif marginal-Efficacité vaccinale
  
############### Calcul des ICs avec l'approche des m-estimateurs ###############

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
    
    # Equations d'estimation : les poids sont estimés à partir d'une régression logistique simple
    
      eq_1 <- (Y == 0)*(V - pscore) # ∂l(β)/β0 = 0
    
      eq_2 <- (Y == 0)*(V - pscore)*C # ∂l(β)/β1 = 0
    
      eq_3 <- (Y*V/plogis(alpha[1] + alpha[2]*C)) - theta[3]
      eq_4 <- (Y*(1 - V)/(1 - plogis(alpha[1] + alpha[2]*C))) - theta[4]
    
      return(c(eq_1, eq_2, eq_3, eq_4))
  }
}

## Estimation avec glm:

  mod.V = glm(V ~ C , data = TNDdat[TNDdat$Y==0,], family = "binomial"); 
  ps = predict(mod.V, newdata = TNDdat, type = "res");
  num = mean((TNDdat$Y*TNDdat$V/ps));
  denom = mean((TNDdat$Y*(1 - TNDdat$V)/(1 - ps))); 

  mestr <- m_estimate(estFUN = geex_ef,                                       
                    data = TNDdat,                                             
                    root_control = setup_root_control(start = c(0, 0, 0.5, 0.5))) # Mêmes valeurs que dans le code que vous m'avez envoyé

  beta_geex <- roots(mestr) # theta = (β0, β1, ψ10)         
  se_geex <- sqrt(diag(vcov(mestr))) # vcov(mestr) : Matrice de variance-covariance

## Comparaison

#c(num, denom) ; beta_geex[3:4]

#RRm <- beta_geex[[3]] / beta_geex[[4]]
#RRm # 0.689678

## Estimation de la variance du RRm
# On suppose que le RRm est distribué selon une loi log-normale,
# theta3 et theta4 sont les valeurs observées des deux paramètres

  var_log_RRm <- (1/theta3^2) * vcov(mestr)[3, 3] + (1/theta4^2) * vcov(mestr)[4, 4]

## Intervalle de confiance pour RRm

  IC_Inf <- theta3/theta4*exp(- 1.96*sqrt(var_log_RRm)) 
  IC_Sup <- theta3/theta4*exp(1.96*sqrt(var_log_RRm)) 

  l <- list(RRm, 1 - RRm, var_log_RRm, IC_Inf, IC_Sup)

  return(l)
   
}