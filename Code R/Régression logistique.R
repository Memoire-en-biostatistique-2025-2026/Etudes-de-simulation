# 1. Calculer les vraies valeurs (fait !)
# 2. Parametres de la simulation
#    - Nombre de replications (1000 ?)
#    - Liste de 1000 germes
# 3. Initialiser des objets pour contenir les resultats
#    - Coefficient de la regression logistique
#    - Erreur-type du coefficient
# 4. Faire une boucle de 1 a nombre de replications
# 5. Generation des donnees TND (n = ?)
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

set.seed(1) # Pour avoir toujours les mêmes germes
nrep <- 10

seeds_list <- sample(1:1000000, size = nrep)

l_vraiRRc <- rep(NA, nrep)
vrai.EV.autres <- rep(NA, nrep)

j = 1;
for (i in seeds_list) {
  
  dat <- datagen.cont(seed = i, popsize = 10000000)
  
  dat0 <- data.frame(C = dat$C, V = 0, Y = dat$I2_0*dat$W2_0*dat$H_0)
  
  dat1 <- data.frame(C = dat$C, V = 1, Y = dat$I2_1*dat$W2_1*dat$H_1) 
  
  dat_complet <- rbind(dat0, dat1)
  
  vraiRRc <- glm(Y ~ V + C, family = binomial(link = "log"), data = dat_complet)
  
  l_vraiRRc[j] <- 1 - exp(coef(vraiRRc)[2])
  
  vraiRRm <- mean(dat1$Y)/mean(dat0$Y)
  
  vrai.EV.autres[j] <- 1 - vraiRRm
  
  j = j + 1
  print(data.frame(Sys.time(), j))
  
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

nrep <- 10

resultats <- data.frame(matrix(ncol = 3, 
                            nrow = nrep))

colnames(resultats) <- c("iteration",
                         "coe_reg", 
                         "err_reg")

# 4. Faire une boucle de 1 a nombre de replications
# 5. Generation des donnees TND (n = 1000)

set.seed(1) # Pour avoir toujours les mêmes germes

l_RRc <- rep(NA, nrep)

seeds_list <- sample(1:1000000, size = nrep)

j <- 1

for (i in seeds_list) {
  
    dat <- datagen(seed = i, ssize = 1000)
    TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
    
# 6. Analyse avec regression logistique
# 7. Enregistrer les resultats
    
    ## TND (infection symptomatique)
    
    fit.TND <- glm(Y ~ V + C,
                       
                       family = binomial(link = "logit"),  
                        # DT: Pour les donnees TND, on utilise logit
                       data = TNDdat) 
    
    resultats.TND <- summary(fit.TND)
    
    ## Numéro de l'itération
    resultats[j, 1] <- j
    
    ## Coefficient de la regression logistique
    resultats[j, 2] <- resultats.TND$coefficients[2]

    ## Erreur-type du coefficient
    resultats[j, 3] <- resultats.TND$coefficients[5]
    
    l_RRc[j] <- 1 - exp(resultats.TND$coefficients[2])
    
    j <- j +1
    
    if(!(j%%10)) print(data.frame(temps = Sys.time(), iter = j))
      # DT : Ajout d'une ligne pour suivre l'avancement

}

# 8. Analyser les resultats

install.packages("Metrics")
library(Metrics)

help(bias)

#    - statistiques descriptives
#    - biais, variance, moyenne de l'erreur-type, couverture des IC

mean(l_RRc)
sd(l_RRc)
  # DT : J'ai remplace RCC par RRc


## Table 01 Résultats de l'étude de simulation (Scénario 01 :
##                                                    P_co-infections = 0.00001)

Tab01 <- data.frame(matrix(ncol = 3, 
                               nrow = 4))

colnames(Tab01) <- c("n",
                      "Methode", 
                      "Regression logistique")

Tab01$n <- c("1000", "-", "-", "-")
Tab01$Methode <- c("Median biais", "MC MSE", "MC SE", "%Cov")

### Calcul du biais

estimates <- 
vrais_valeurs <-   

biais <- 
  
mean_bias <- mean(bias, na.rm = TRUE)
median_bias <- median(bias, na.rm = TRUE)
rmse <- sqrt(mean(bias^2, na.rm = TRUE))
mc_standard_error <- sd(estimates, na.rm = TRUE) / sqrt(length(estimates))
