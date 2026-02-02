
source("Code R/Génération des données.R")
source("Code R/Fonc01_RegLog.R")
source("Code R/Fonc02_IPW.R")
source("Code R/Fonc03_TNDDR.R")

# Chargements des librairies nécessaires

library(simhelpers)
library(dplyr)
library(tibble)
library(knitr)
library(dplyr)
library(kableExtra)

################################################################################
# Scénario 02 (co-infection ~ 10% dans l'échantillon): Sous-scénario 01 : CV = 33%, 
#                                                                         I1_prev = 15% et 
#                                                                         I2_rev = 50%
################################################################################

# Calcul des vraies valeurs 

set.seed(1) # Pour avoir toujours les mêmes germes

nsim <- 10 # Nombre de réplications

seeds_list <- sample(1:1000000, size = nsim)

l_vraiRRc <- rep(NA, nsim)
l_vraiRRm <- rep(NA, nsim)

for (i in 1:nsim) {
  
  dat <- datagen.cont(seed = seeds_list[i], popsize = 1000000,
                      co_inf_para1 = 2, co_inf_para2 = -2)
  
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
sd(l_vraiRRc)


l_vraiRRm

mean(l_vraiRRm)
sd(l_vraiRRm)

################################################################################
# Initialiser des objets pour contenir les resultats

Tab01 <- data.frame(n = c("1000", "-", "-", "-"))

################################################################################
########################## Analyse des résultats ###############################

nsim <- 100

# Listes des valeurs à tester (à discuter)

CV_liste <- c(0.33, 0.5, 0.7)
I1_prev_liste <- c(0.15, 0.2, 0.5)
I2_prev_liste <- c(0.5, 0.2, 0.15)

# Initialiser des objets pour contenir les resultats

resultats <- data.frame(matrix(ncol = 4, 
                               nrow = nsim * 3))

colnames(resultats) <- c("coe_reg", 
                         "err_reg",
                         "RRc", # Risque relatif conditionnel
                         "est_VE")

resultats2 <- data.frame(matrix(ncol = 5, 
                                nrow = nsim * 3))

colnames(resultats2) <- c("RRm",# Risque relatif marginal
                          "VE",
                          "var_log_RRm",# Variance du log du risque relatif marginal
                          "IC_inf", # Borne inférieure de l'intervalle de confiance
                          "IC_sup") # Sa borne supérieure

resultats3 <- data.frame(matrix(ncol = 36, 
                                nrow = nsim * 3))

colnames(resultats3) <- c("RRm_RF", "VE_RF", "var_log_RRm-RF", "IC_inf1-RF", "IC_sup1-RF", "IC_inf2-RF", "IC_sup2-RF", 
                          "IC_inf3-RF", "IC_sup3-RF",# Forêt aléatoire
                          
                          "RRm_Lasso", "VE_Lasso","var_log_RRm-Lasso", "IC_inf1-Lasso", "IC_sup1-Lasso", 
                          "IC_inf2-Lasso", "IC_sup2-Lasso", 
                          "IC_inf3-Lasso", "IC_sup3-Lasso",# Régression Lasso
                          
                          "RRm_Mars", "VE_Mars","var_log_RRm-Mars", "IC_inf1-Mars", "IC_sup1-Mars", 
                          "IC_inf2-Mars", "IC_sup2-Mars", 
                          "IC_inf3-Mars", "IC_sup3-Mars",# Régression avec splines
                          
                          "RRm_RN", "VE_RN","var_log_RRm-RN", "IC_inf1-RN", "IC_sup1-RN", 
                          "IC_inf2-RN", "IC_sup2-RN", 
                          "IC_inf3-RN", "IC_sup3-RN"# Réseaux de neurones
                          
                          
                          # Risque relatid marginal   
                          # Variance du log du risque relatif marginal
                          # Borne inférieure du premier intervalle de confiance
                          # Sa borne supérieure
                          # Borne inférieure du deuxième intervalle de confiance
                          # Sa borne supérieure
                          # Borne inférieure du troisième intervalle de confiance
                          # Sa borne supérieure
                          
)

methode <- list(RandomForest, Lasso, Mars, RN)

############################### Estimation #####################################

d <- 1
f <- nsim

for(j in CV_liste) {
  
  for (i in d:f) {
    
  
    dat <- datagen(seed = seeds_list[i], ssize = 1000, co_inf_para1 = 2, co_inf_para2 = -2, CV = j, popsize = 1*10**6)
  
    resultats[i,] <- RegLog(dat) # Régression logistique
    resultats2[i,] <- IPW(dat)   # IPW
  
    l <- list()
  
    for(k in methode) {
      
      m <- TNDDR(dat, k) # Liste des résultats pour la méthode k
      l <- append(l, m) # Combiner les résultats des différentes méthodes
    
  }
  
    resultats3[i,] <- l
  
  
  # DT : Ajout d'une ligne pour suivre l'avancement
  
    if(!(i%%10)) print(data.frame(temps = Sys.time(), iter = i))
  
  }
  
  d <- d + nsim
  f <- f + nsim
  
}


########################## Affichage des résultats #############################


