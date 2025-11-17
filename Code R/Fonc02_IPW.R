
# Chargement des librairies nécessaires

library(simhelpers)
library(dplyr)
library(tibble)
library(knitr)
library(dplyr)
library(kableExtra)
library(geepack)

# Initialiser des objets pour contenir les resultats

resultats2 <- data.frame(matrix(ncol = 4, 
                                nrow = nsim))

colnames(resultats2) <- c("iteration",
                          "germe",
                          "RRm_IPW",
                          "VE") # Risque relatif marginal

seeds_list <- sample(1:1000000, size = nsim)

for (i in 1:nsim) {
  
  dat <- datagen(seed = seeds_list[i], ssize = 1000)
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

##    - statistiques descriptives

summary(resultats2$RRm)
mean(resultats2$RRm)
sd(resultats2$RRm)

##    - biais, variance, moyenne de l'erreur-type, couverture des IC

resultats2$vrai_param <- rep(0.6326075, nsim) # vraie valeur du paramètre

Tab01 <- data.frame(matrix(ncol = 3, 
                           nrow = 4))

colnames(Tab01) <- c("n",
                     "Methode", 
                     "IPW")

Tab01$n <- c("1000", "-", "-", "-")
Tab01$Methode <- c("MCSE_bias", "MCSE_var", "MCSE_mse", "%Cov")

help("calc_absolute") # calculer les différentes mesures de performance

### MCSE_biais

MCSE_biais <- calc_absolute(resultats2, RRm, vrai_param, criteria = "bias")
Tab01$IPW[1] <- MCSE_biais[3]

### MCSE_var

MCSE_var <- calc_absolute(resultats2, RRm, vrai_param, criteria = "var")
Tab01$IPW[2] <- MCSE_var[3]


### MCSE_MSE 

MCSE_MSE <- calc_absolute(resultats2, RRm, vrai_param, criteria = "mse")
Tab01$IPW[3] <- MCSE_MSE[3]

### Coverage

################################################################################
############### Calcul des ICs avec l'approche des m-estimateurs ###############

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

help("m_estimate")

mestr <- m_estimate(estFUN = geex_ef,                                       
                    data = TNDdat,                                             
                    root_control = setup_root_control(start = c(0, 0, 0.5, 0.5))) # Mêmes valeurs que dans le code que vous m'avez envoyé

beta_geex <- roots(mestr) # theta = (β0, β1, ψ10)         
se_geex <- sqrt(diag(vcov(mestr))) # vcov(mestr) : Matrice de variance-covariance

## Comparaison

c(num, denom) ; beta_geex[3:4]

RRm <- beta_geex[[3]] / beta_geex[[4]]
RRm # 0.6245755

## Estimation de la variance du mRR
# On suppose que le RRm est distribué selon une loi log-normale,
# theta3 et theta4 sont les valeurs observées des deux paramètres

var_log_RRm <- (1/theta3^2) * vcov(mestr)[3, 3] + (1/theta4^2) * vcov(mestr)[4, 4]

## Intervalle de confiance pour RRm

IC_Inf <- theta3/theta4*exp(- 1.96*sqrt(var_log_RRm)) 
IC_Sup <- theta3/theta4*exp(1.96*sqrt(var_log_RRm)) 

kable(Tab01)