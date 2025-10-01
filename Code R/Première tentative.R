# Ce fichier va contenir ma première tentative de génération de données basée sur le DAG que vous m'avez envoyé

# Génération des données

## Chargement des librairies nécessaires

library(simDAG)
library(ggplot2)
library(ggforce) # Pour afficher le DAG pour la deuxième méthode de génération

datagen <- function(seed = sample(1:1000000, size = 1), ssize = 5000, 
                    
                    popsize = 15000000) {
  
          
          # Génération du facteur de confusion continu C 
          
          C <- runif(n = popsize, 0.1, 3); # On commence par les noeuds racines de notre DAG
          
          # Génération du statut vaccinal V : V ~ Bernoulli(logit(a0 + a1*C))
          
            
            V <- rbinom(n = popsize, size = 1, prob = plogis(-12.25 + 0.3*C))
          
          # Génération des infections I1 : I1 ~ Bernoulli(logit(a0 + a1*C)) et 
          #                           I2 : I2 ~ Bernoulli(logit(a0 + a1*C + a2*V))
          

          I1 <- rbinom(n = popsize, size = 1, prob = plogis(-11.51292 + 0.35*C))  
          
          I2 <- rbinom(n = popsize, size = 1, prob = plogis(-11.51292 + 0.15*C - 0.1*V))  
          
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
################################################################################

dag <- empty_dag() +
        
       node("C", type = "runif", 0.1, 3) +
       
       node("V", type = "binomial", formula = ~ -12.25 + 0.3*C) +
       
       node("I1", type = "binomial", formula = ~ -11.51292 + 0.35*C) 

  
set.seed(10)

sim_dat <- sim_from_dag(dag, n_sim = 100000)

plot(dag)

  
  