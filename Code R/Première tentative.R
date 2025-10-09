# Ce fichier va contenir ma première tentative de génération de données basée sur le DAG que vous m'avez envoyé
# Les différentes valeurs des paramètres utilisées sont issues de la lescture des différents articles 

# Génération des données

## Chargement des librairies nécessaires

library(simDAG)
library(ggplot2)
library(ggforce) # Pour afficher le DAG pour la deuxième méthode de génération

datagen <- function(seed = sample(1:1000000, size = 1), ssize = 5000, 
                         
                         popsize = 150000) {
  
  
  # Génération du facteur de confusion continu C 
  
  C <- runif(n = popsize, 0.1, 3); # On commence par les noeuds racines de notre DAG
  
  # Génération des infections I1 : I1 ~ Bernoulli(logit(a0 + a1*C)) et 
  #                           I2 : I2 ~ Bernoulli(logit(a0 + a1*C + a2*V))
  
  
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(-2.5 + 0.35*C))
  
  I2 <- rbinom(n = popsize, size = 1, prob = plogis(0 + 0.15*C - 0.1*V))
  
  
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
  
  # Génération de W
  
  W <- pmax(W1, W2)
  
  # Génération de l'hospitalization 
  
  H = rep(0, popsize) # H = 0 si W = 0, donc:

  H[W == 1] <- rbinom(prob = plogis(-1.5 + 0.5*C[W == 1]),
                          size = 1, n = sum(W == 1))
  
  # Le devis test-négatif ne conserve que les personnes testées : H = 1
  
  R <- sample(which(H == 1), ssize, replace = TRUE) 
  
  if (return_full == FALSE) {
    
    dat <- as.data.frame(cbind(Y = I2, V = V, C = C)[R, ])
    
  } else {
    
    dat <- as.data.frame(cbind(Infec_RSV = I2, Infec = I1, H = H, W = W, V = V,
                               C = C)) # Virus respiratoire syncytial RSV
    
  }
  
  return(dat)
  
}  

dat <- datagen(ssize = 5000)

################################################################################ 

datagen.cont <- function(seed = sample(1:1000000, size = 1), ssize = 5000, 
                    
                         popsize = 150000) {
  
  
  # Génération du facteur de confusion continu C 
  
  C <- runif(n = popsize, 0.1, 3); # On commence par les noeuds racines de notre DAG
  
  # Génération des infections I1 : I1 ~ Bernoulli(logit(a0 + a1*C)) et 
  #                           I2 : I2 ~ Bernoulli(logit(a0 + a1*C + a2*V))
  
  
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(-2.5 + 0.35*C))
  
  I2_0 <- rbinom(n = popsize, size = 1, prob = plogis(0 + 0.15*C - 0.1*0))
  I2_1 <- rbinom(n = popsize, size = 1, prob = plogis(0 + 0.15*C - 0.1*1))
  
  
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
}

################################################################################

datagen_con_2 <- function(seed = sample(1:1000000, size = 1), ssize = 5000, 
                    
                    popsize = 150000, cfV0 = FALSE, cfV1 = FALSE) {
  
          
          # Génération du facteur de confusion continu C 
          
          C <- runif(n = popsize, 0.1, 3); # On commence par les noeuds racines de notre DAG
          
          # Génération du statut vaccinal V : V ~ Bernoulli(logit(a0 + a1*C))
          
          if (cfV0 == TRUE) {
            
            V <- rep(0, popsize)
            
          } else if (cfV1 == TRUE) {
            
            V <- rep(1, popsize)
            
          } else {
            
          V <- rbinom(n = popsize, size = 1, prob = plogis(0.5 + 0.3*C))

          }
          
          # Génération des infections I1 : I1 ~ Bernoulli(logit(a0 + a1*C)) et 
          #                           I2 : I2 ~ Bernoulli(logit(a0 + a1*C + a2*V))
          

          I1 <- rbinom(n = popsize, size = 1, prob = plogis(-2.5 + 0.35*C))

          
          I2 <- rbinom(n = popsize, size = 1, prob = plogis(0 + 0.15*C - 0.1*V))

          
          # Génération des symptomes W1: W1~Bernoulli(logit(a0 + a1*C[I1 = 1]))  et 
          #                          W2: W2~Bernoulli(logit(a0 + a1*C[I2 == 1] + a2*V[I2 == 1]))
          
          W1 <- rep(0, popsize) # W1 = 0 si I1 = 0, donc:
          
          W1[I1 == 1] <- rbinom(
            
            n = sum(I1 == 1),
            size = 1,
            prob = plogis(-0.5 + 0.5 * C[I1 == 1])
            
          )
          
          W2 <- rep(0, popsize) # W2 = 0 si I2 = 0, donc:
          
          W2[I2 == 1] <- rbinom(
            
            n = sum(I2 == 1),
            size = 1,
            prob = plogis(-3.75 + 2*C[I2 == 1] - 0.91*V[I2 == 1])
            
          )
          
          # Génération de W
          
          W <- pmax(W1, W2)
          
          # Génération de l'hospitalization 
          
          H = rep(0, popsize) # H = 0 si W = 0, donc:
          
          H[W == 1] <- rbinom(prob = plogis(-1.5 + 0.5*C[W == 1]),
                              size = 1, n = sum(W == 1))
   

}
