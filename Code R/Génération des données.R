
## Chargement des librairies nécessaires

library("stats")
library("rje")


datagen <- function(seed = sample(1:1000000, size = 1), ssize = 5000, 
                    
                    popsize = 150000, co_inf_para1 = 0, co_inf_para2 = 0,
                    
                    CV = 0.33, # Couverture vaccinale
                    I1_prev = 0.15, # Prévalence de I1
                    I2_prev = 0.5, # Prévalence de I2
                    W1_prev = 0.05,# Prévalence de W1
                    W2_prev = 0.02, # Prévalence de W2
                                   # Valeurs par défaut pour le scénario de base
                    
                    return_full = FALSE) { # Pour choisir entre population et échantillon
  
  
  # Génération du facteur de confusion continu C 
  
  C <- runif(n = popsize, 0.1, 3) # On commence par les noeuds racines de notre DAG
  
  # Génération du statut vaccinal V : V ~ Bernoulli(logit(a0 + a1*C))
  
  a <- 1/0.87
  # Définir la fonction en b0 (en calculant l'intégrale) pour V
  
  f_V <- function (b0) {
    
    a*(log(1 + exp(b0 + 3*0.3)) - log(1 + exp(b0 + 0.1*0.3))) - CV
      
       # Contrainte pour la probabilité marginale
    
  }
  
  b0 <- uniroot(f_V, c(-10, 10))$root
  
  V <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.3*C)) # Couverture vaccinale (hyperparamètre)
  # Voir fichier Anciens fichiers\Génération des variables
  # pour la valeur de b0 
  
  # Génération des infections I1 : I1 ~ Bernoulli(logit(b0 + b1*C + b2*C2)) et 
  #                           I2 : I2 ~ Bernoulli(logit(b0 + b1*C + b2*V + b3*C2))
  
  C2 <- 2*rbinom(n = popsize, 1, 0.5) - 1
  
  ## Calcul de la valeur de b0 pour I1 avec b1 et b2 connus
  
  # On pose:
  
  a <- 0.5/2.9

  # P(I1 = 1) - 0.15
  
  f_I1 <- function (b0) {
    
    (a/0.35)*(log(1 + exp(b0 + 3*0.35 - co_inf_para1)) + log(1 + exp(b0 + 3*0.35 + co_inf_para1)) - 
              
            log(1 + exp(b0 + 0.1*0.35 - co_inf_para1)) - log(1 + exp(b0 + 0.1*0.35 + co_inf_para1))) - 
      
      I1_prev # Contrainte pour la probabilité marginale
    
  }
  
  b0 <- uniroot(f_I1, c(-10, 10))$root
  
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.35*C + co_inf_para1*C2)) # Pour une prévalence ~ 
  
  ## Calcul de la valeur de b0 pour I2 avec b1, b2 et b3 connus
  
  P_I2 <- function(c,b0) { # Probabilité conjointe pour I2
    
    a*(expit(b0 + 0.15*c - 0.1 - co_inf_para2)*expit(0.5 + 0.3*c) + 
         
         expit(b0 + 0.15*c - co_inf_para2)*(1 - expit(0.5 + 0.3*c)) + 
         
         expit(b0 + 0.15*c - 0.1 + co_inf_para2)*expit(0.5 + 0.3*c) + 
         
         expit(b0 + 0.15*c + co_inf_para2)*(1 - expit(0.5 + 0.3*c))
    )
    
  }
  
  # P(I2 = 2) - 0.50
  
  f_I2 <- function(b0){
    
    integrate(P_I2, 0.1, 3, b0)$value - I2_prev # Contrainte pour la probabilité marginale
    
  }
  
  b0 <- uniroot(f_I2, c(-10, 10))$root
  
  I2 <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.15*C - 0.1*V + co_inf_para2*C2)) # Pour une prévalence ~ 
  
  # Génération des symptomes W1: W1~Bernoulli(logit(a0 + a1*C[I1 = 1]))  et 
  #                          W2: W2~Bernoulli(logit(a0 + a1*C[I2 == 1] + a2*V[I2 == 1]))
  
  # Définir la fonction en b0 (en calculant l'intégrale) pour W1
  
  a <- 1/2.9
  b1 <- 0.5
  
  f_W1 <- function (b0) {
    
    (a/b1)*(log(1 + exp(b0 + 3*b1 ))  - log(1 + exp(b0 + 0.1*b1))) - 
      
      W1_prev # Contrainte pour la probabilité marginale
    
  }
  
  b0 <- uniroot(f_W1, c(-10, 10))$root 
  
  W1 <- rep(0, popsize) # W1 = 0 si I1 = 0, donc:
  
  W1[I1 == 1] <- rbinom(
    
    n = sum(I1 == 1),
    size = 1,
    prob = plogis(b0 + 0.5 * C[I1 == 1])
    
  ) # W1 = 1 ~ 4%-5%
  
  # Définir la fonction en b0 (en calculant l'intégrale) pour W2
  
  # On pose:
  
  a <- 0.5/2.9
  a0 <- -1.18
  a1 <- 0.3
  b1 <- 2
  b2 <- -0.91
  
  P_W2 <- function(c,b0) { # Probabilité conjointe pour W2
    
    a*(expit(b0 + b1*c + b2)*expit(a0 + a1*c) + 
         
         expit(b0 + b1*c)*expit(a0 + a1*c) + 
         
         expit(b0 + b1*c + b2)*expit(a0 + a1*c) + 
         
         expit(b0 + b1*c )*expit(a0 + a1*c)
    )
    
  }
  
  help(integrate)
  # ...	
  # additional arguments to be passed to f : on utilise cet argument pour passer 
  # la valeur b0 à la fonction P_I2 dans ce cas.
  
  f_W2 <- function(b0){
    
    integrate(P_W2, 0.1, 3, b0)$value - W2_prev # Contrainte pour la probabilité marginale
    
  }
  
  b0 <- uniroot(f_W2, c(-10, 10))$root 
   
  W2 <- rep(0, popsize) # W2_0 = 0 si I2_0 = 0, donc:
  
  W2[I2 == 1] <- rbinom(
    
    n = sum(I2 == 1),
    size = 1,
    prob = plogis(b0 + 2*C[I2 == 1] - 0.91*V[I2 == 1])
    
  ) # W2 = 1 ~ 3%
  
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

# dat <- datagen(ssize = 5000, co_inf_para1 = 8, co_inf_para2 = -8)
dat_full <- datagen(return_full = TRUE)

sum(dat_full$V == 1)/nrow(dat_full) * 100 # Couverture vaccinale

sum(dat_full$Infec == 1)/nrow(dat_full) * 100
sum(dat_full$W1 == 1)/sum(dat_full$Infec == 1) * 100 # Prévalence symptomatique pour I1

sum(dat_full$Infec_RSV == 1)/nrow(dat_full) * 100
sum(dat_full$W2 == 1)/sum(dat_full$Infec_RSV == 1) * 100 # Prévalence symptomatique pour I2

################################################################################ 

# Génération des scénarios contrefactuels

datagen.cont <- function(seed = sample(1:1000000, size = 1), popsize = 150000,
                         
                         co_inf_para1 = 0, co_inf_para2 = 0,
                         
                         CV = 0.33, # Couverture vaccinale
                         I1_prev = 0.15, # Prévalence de I1
                         I2_prev = 0.5,  # Prévalence de I2
                         W1_prev = 0.05, # Prévalence de W1
                         W2_prev = 0.02 # Prévalence de W2
                                         # Valeurs par défaut pour le scénario de base
                         ) {
  
  
  # Génération du facteur de confusion continu C 
  
  C <- runif(n = popsize, 0.1, 3); # On commence par les noeuds racines de notre DAG
  
  # Génération du statut vaccinal V : V ~ Bernoulli(logit(b0 + b1*C))
  
  a <- 1/0.87
  # Définir la fonction en b0 (en calculant l'intégrale) pour V
  
  f_V <- function (b0) {
    
    a*(log(1 + exp(b0 + 3*0.3)) - log(1 + exp(b0 + 0.1*0.3))) - 
      
      CV # Contrainte pour la probabilité marginale
    
  }
  
  b0 <- uniroot(f_V, c(-10, 10))$root
  
  V <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.3*C)) # Couverture vaccinale (hyperparamètre)
  # Voir fichier Anciens fichiers\Génération des variables
  # pour la valeur de b0 
  
  # Génération des infections I1 : I1 ~ Bernoulli(logit(a0 + a1*C)) et 
  #                           I2 : I2 ~ Bernoulli(logit(a0 + a1*C + a2*V))
  
  
  C2 <- 2*rbinom(n = popsize, 1, 0.5) - 1
  
  ## Calcul de la valeur de b0 pour I1 avec b1 et b2 connus
  
  # On pose:
  
  a <- 0.5/2.9
  
  # P(I1 = 1) - 0.15
  
  f_I1 <- function (b0) {
    
    (a/0.35)*(log(1 + exp(b0 + 3*0.35 - co_inf_para1)) + log(1 + exp(b0 + 3*0.35 + co_inf_para1)) - 
                
              log(1 + exp(b0 + 0.1*0.35 - co_inf_para1)) - log(1 + exp(b0 + 0.1*0.35 + co_inf_para1))) - 
      
      I1_prev # Contrainte pour la probabilité marginale
    
  }
  
  b0 <- uniroot(f_I1, c(-10, 10))$root
  
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.35*C + co_inf_para1*C2))
  
  ## Calcul de la valeur de b0 pour I2 avec b1, b2 et b3 connus
  
  P_I2 <- function(c,b0) { # Probabilité conjointe pour I2
    
    a*(expit(b0 + 0.15*c - 0.1 - co_inf_para2)*expit(0.5 + 0.3*c) + 
         
         expit(b0 + 0.15*c - co_inf_para2)*(1 - expit(0.5 + 0.3*c)) + 
         
         expit(b0 + 0.15*c - 0.1 + co_inf_para2)*expit(0.5 + 0.3*c) + 
         
         expit(b0 + 0.15*c + co_inf_para2)*(1 - expit(0.5 + 0.3*c))
    )
    
  }
  
  # P(I2 = 2) - 0.50
  
  f_I2 <- function(b0){
    
    integrate(P_I2, 0.1, 3, b0)$value - I2_prev # Contrainte pour la probabilité marginale
    
  }
  
  b0 <- uniroot(f_I2, c(-10, 10))$root
  
  I2_0 <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.15*C - 0.1*0 + co_inf_para2*C2))
  I2_1 <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.15*C - 0.1*1 + co_inf_para2*C2))
  
  # Génération des symptomes W1: W1~Bernoulli(logit(a0 + a1*C[I1 = 1]))  et 
  #                          W2: W2~Bernoulli(logit(a0 + a1*C[I2 == 1] + a2*V[I2 == 1]))
  
  # Définir la fonction en b0 (en calculant l'intégrale) pour W1
  
  a <- 1/2.9
  b1 <- 0.5
  
  f_W1 <- function (b0) {
    
    (a/b1)*(log(1 + exp(b0 + 3*b1 ))  - log(1 + exp(b0 + 0.1*b1))) - 
      
      W1_prev # Contrainte pour la probabilité marginale
    
  }
  
  b0 <- uniroot(f_W1, c(-10, 10))$root
  W1 <- rep(0, popsize) # W1 = 0 si I1 = 0, donc:
  
  W1[I1 == 1] <- rbinom(
    
    n = sum(I1 == 1),
    size = 1,
    prob = plogis(b0 + 0.5 * C[I1 == 1])
    
  )
  
  # Définir la fonction en b0 (en calculant l'intégrale) pour W2
  
  # On pose:
  
  a <- 0.5/2.9
  a0 <- -1.18
  a1 <- 0.3
  b1 <- 2
  b2 <- -0.91
  
  P_W2 <- function(c,b0) { # Probabilité conjointe pour W2
    
    a*(expit(b0 + b1*c + b2)*expit(a0 + a1*c) + 
         
         expit(b0 + b1*c)*expit(a0 + a1*c) + 
         
         expit(b0 + b1*c + b2)*expit(a0 + a1*c) + 
         
         expit(b0 + b1*c )*expit(a0 + a1*c)
    )
    
  }
  
  help(integrate)
  # ...	
  # additional arguments to be passed to f : on utilise cet argument pour passer 
  # la valeur b0 à la fonction P_I2 dans ce cas.
  
  f_W2 <- function(b0){
    
    integrate(P_W2, 0.1, 3, b0)$value - W2_prev # Contrainte pour la probabilité marginale
    
  }
  
  b0 <- uniroot(f_W2, c(-10, 10))$root
  
  W2_0 <- rep(0, popsize) # W2_0 = 0 si I2_0 = 0, donc:
  W2_1 <- rep(0, popsize) # W2_1 = 0 si I2_1 = 0, donc:
  
  W2_0[I2_0 == 1] <- rbinom(
    
    n = sum(I2_0 == 1),
    size = 1,
    prob = plogis(b0 + 2*C[I2_0 == 1] - 0.91*0) # V[I2_0 == 1]
    
  )
  
  W2_1[I2_1 == 1] <- rbinom(
    
    n = sum(I2_1 == 1), 
    size = 1,
    prob = plogis(b0 + 2*C[I2_1 == 1] - 0.91*1) # V[I2_1 == 1]
    
  )
  
  # Génération de W
  
  W_0 <- pmax(W1, W2_0)
  W_1 <- pmax(W1, W2_1)
  
  
  # Génération de l'hospitalization 
  
  H_0 = rep(0, popsize) # H_0 = 0 si W_0 = 0, donc:
  H_1 = rep(0, popsize) # H_1 = 0 si W_1 = 0, donc:
  
  H_0[W_0 == 1] <- rbinom(prob = plogis(-1.5 + 0.5*C[W_0 == 1]),
                          size = 1, n = sum(W_0 == 1))
  
  H_1[W_1 == 1] <- rbinom(prob = plogis(-1.5 + 0.5*C[W_1 == 1]),
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

# nsim <- 10 # Fixer le nombre de réplications
