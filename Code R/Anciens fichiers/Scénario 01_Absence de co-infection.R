# Scénario 01 ~ Absence de co-infection

# Un scénario où on s'attend à ce que les deux méthodes d'estimation fonctionnent
# correctement pour voir s'il n'ya pas d'erreur dans le code.

# En fixant a0 pour I1 et I2 à une valeur qui permet d'obtenir de faibles taux de 
# co-infection, de sorte que l'hypothèse d'échangeabilité de contrôle se vérifie approximativement

datagen <- function(seed = sample(1:1000000, size = 1), ssize = 5000, 
                    
                    popsize = 150000, co_inf_para1 = -11, co_inf_para2 = 11,
                    
                    return_full = FALSE) { # Pour choisir entre population et échantillon
  
  
  # Génération du facteur de confusion continu C 
  
  C <- runif(n = popsize, 0.1, 3) # On commence par les noeuds racines de notre DAG
  
  # Génération du statut vaccinal V : V ~ Bernoulli(logit(a0 + a1*C))
  
  V <- rbinom(n = popsize, size = 1, prob = plogis(0.5 + 0.3*C))
  
  # Génération des infections I1 : I1 ~ Bernoulli(logit(a0 + a1*C)) et 
  #                           I2 : I2 ~ Bernoulli(logit(a0 + a1*C + a2*V))
  
  C2 <- 2*rbinom(n = popsize, 1, 0.5) - 1
  
  p1_temp = plogis(0.35*C + co_inf_para1*C2);
  p1 = p1_temp/mean(p1_temp)*0.15;
  
  I1 <- rbinom(n = popsize, size = 1, prob = p1)
  
  p2_temp = plogis(0 + 0.15*C - 0.1*V + co_inf_para2*C2)
  p2 = p2_temp/mean(p2_temp)*0.50;
  
  I2 <- rbinom(n = popsize, size = 1, prob = p2)
  
  # Génération des symptomes W1: W1~Bernoulli(logit(a0 + a1*C[I1 = 1]))  et 
  #                          W2: W2~Bernoulli(logit(a0 + a1*C[I2 == 1] + a2*V[I2 == 1]))
  
  W1 <- rep(0, popsize) # W1 = 0 si I1 = 0, donc:
  
  W1[I1 == 1] <- rbinom(
    
    n = sum(I1 == 1),
    size = 1,
    prob = plogis(-0.5 + 0.5 * C[I1 == 1])
    
  )
  
  W2 <- rep(0, popsize) # W2_0 = 0 si I2_0 = 0, donc:
  
  W2[I2 == 1] <- rbinom(
    
    n = sum(I2 == 1),
    size = 1,
    prob = plogis(-3.75 + 2*C[I2 == 1] - 0.91*V[I2 == 1])
    
  )
  
  co_W <- sum(W1 == 1 & W2 == 1, na.rm = TRUE)
  per_co_W <- co_W / popsize * 100
  
  # Génération de W
  
  W <- pmax(W1, W2)
  
  # Génération de l'hospitalization 
  
  H = rep(0, popsize) # H = 0 si W = 0, donc:
  
  H[W == 1] <- rbinom(prob = plogis(-1.5 + 0.5*C[W == 1]),
                      size = 1, n = sum(W == 1))
  
  # Le devis test-négatif ne conserve que les personnes testées : H = 1
  
  R <- sample(which(H == 1), ssize, replace = TRUE) # Échantillon aléatoire parmi les personnes hospitalisées
  
  if (return_full == FALSE) {
    
    dat <- as.data.frame(cbind(Infec_RSV = I2, Infec = I1, H = H, W1 = W1, W2 = W2,
                               
                               W = W, V = V, C = C)[R, ]) # Virus respiratoire syncytial RSV
    
    # Calcul du pourcentage de co-infection (symptomatique) dans l'échantillon
    
    co_inf <- sum(dat$Infec == 1 & dat$Infec_RSV == 1)
    
    per_co_inf <- co_inf / ssize * 100
    
    print(paste("Le pourcentage de co_infection dans l'échantillon est :", per_co_inf))
    
  } else { # Données pour la population totale (return_full == TRUE)
    
    dat <- as.data.frame(cbind(Infec_RSV = I2, Infec = I1, H = H, W1 = W1, W2 = W2,
                               
                               W = W, V = V, C = C))
    
    # Calcul du pourcentage de co-infection symptomatique dans la population
    
    co_inf_1 <- sum(dat$Infec == 1 & dat$Infec_RSV == 1 & dat$W == 1)
    
    per_co_inf_1 <- co_inf_1 / popsize * 100
    
    print((paste("Le pourcentage de co_infection symptomatique dans la population est :", per_co_inf_1)))
    
    # Calcul du pourcentage de co-infection asymptomatique dans la population
    
    co_inf_0 <- sum(dat$Infec == 1 & dat$Infec_RSV == 1 & dat$W == 0)
    
    per_co_inf_0 <- co_inf_0 / popsize * 100
    
    print((paste("Le pourcentage de co_infection asymptomatique dans la population est :", per_co_inf_0)))
    
  }
  
  return(dat)
  
}  

dat <- datagen(ssize = 5000)
#  "Le pourcentage de co_infection dans l'échantillon est : 0"

dat_full <- datagen(return_full = TRUE)
# "Le pourcentage de co_infection symptomatique dans la population est : 0"
# "Le pourcentage de co_infection asymptomatique dans la population est : 0"

################################################################################

datagen.cont <- function(seed = sample(1:1000000, size = 1), popsize = 150000,
                         
                         co_inf_para = 0) {
  
  
  # Génération du facteur de confusion continu C 
  
  C <- runif(n = popsize, 0.1, 3); # On commence par les noeuds racines de notre DAG
  
  # Génération des infections I1 : I1 ~ Bernoulli(logit(a0 + a1*C)) et 
  #                           I2 : I2 ~ Bernoulli(logit(a0 + a1*C + a2*V))
  
  
  C2 <- runif(n = popsize, 0,1)
  
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(-7 + 0.35*C + co_inf_para*C2))
  
  I2_0 <- rbinom(n = popsize, size = 1, prob = plogis(-7 + 0.15*C + co_inf_para*C2 - 0.1*0))
  I2_1 <- rbinom(n = popsize, size = 1, prob = plogis(-7 + 0.15*C + co_inf_para*C2 - 0.1*1))
  
  # Génération des symptomes W1: W1~Bernoulli(logit(a0 + a1*C[I1 = 1]))  et 
  #                          W2: W2~Bernoulli(logit(a0 + a1*C[I2 == 1] + a2*V[I2 == 1]))
  
  W1 <- rep(0, popsize) # W1 = 0 si I1 = 0, donc:
  
  W1[I1 == 1] <- rbinom(
    
    n = sum(I1 == 1),
    size = 1,
    prob = plogis(-0.5 + 0.5 * C[I1 == 1])
    
  )
  
  W2_0 <- rep(0, popsize) # W2_0 = 0 si I2_0 = 0, donc:
  W2_1 <- rep(0, popsize) # W2_1 = 0 si I2_1 = 0, donc:
  
  W2_0[I2_0 == 1] <- rbinom(
    
    n = sum(I2_0 == 1),
    size = 1,
    prob = plogis(-3.75 + 2*C[I2_0 == 1] - 0.91*0) # V[I2_0 == 1]
    
  )
  
  W2_1[I2_1 == 1] <- rbinom(
    
    n = sum(I2_1 == 1), 
    size = 1,
    prob = plogis(-3.75 + 2*C[I2_1 == 1] - 0.91*1) # V[I2_1 == 1]
    
  )
  
  # Génération de W
  
  W_0 <- pmax(W1, W2_0)
  W_1 <- pmax(W1, W2_1)
  
  
  # Génération de l'hospitalization 
  
  H_0 = rep(0, popsize) # H_0 = 0 si W_0 = 0, donc:
  H_1 = rep(0, popsize) # H_1 = 0 si W_1 = 0, donc:
  
  H_0[W_0 == 1] <- rbinom(prob = plogis(-1.5 + 0.5*C[W_0 == 1]),
                          size = 1, n = sum(W_0 == 1))
  
  H_1[W_1 == 1] <- rbinom(prob = plogis(-1.5 + 0.5*C[W_0 == 1]),
                          size = 1, n = sum(W_1 == 1))
  
  dat <- data.frame(C, I1, I2_0, I2_1, W1, W2_0, W2_1, W_0, W_1, H_0, H_1)
  
  # Calcul du pourcentage de co-infection
  
  co_inf_0 <- sum(dat$I1 == 1 & dat$I2_0 == 1)
  
  per_co_inf_0 <- co_inf_0 / popsize * 100
  
  co_inf_1 <- sum(dat$I1 == 1 & dat$I2_1 == 1)
  
  per_co_inf_1 <- co_inf_1 / popsize * 100
  
  print(paste(per_co_inf_0, per_co_inf_1))
  
  return(dat)
}

################################################################################

### Calcul des vraies valeurs ### popsize = 1000000!

## Regression logistique ##

l_vraiVE
# 0.4129067 0.4103635 0.5108003 0.3956070 0.4068455 0.3146842 0.4200826 0.3269186 0.3796457 0.4786688

mean(l_vraiVE)
# 0.4056523

sd(l_vraiVE)
# 0.05955708

# |n    |Methode   |Regression_logistique |
# |:----|:---------|:---------------------|
# |1000 |MCSE_bias |0.0632358             |
# |-    |MCSE_var  |0.1626116             |
# |-    |MCSE_mse  |0.1818748             |
# |-    |%Cov      |0.38                  |
