# Ce fichier va contenir ma première tentative de génération de données basée sur le DAG que vous m'avez envoyé

# Génération des données

## Chargement des librairies nécessaires

library(simDAG)
library(ggplot2)
library(ggforce) # Pour afficher le DAG

datagen <- function(seed = sample(1:1000000, size = 1), ssize = 5000, 
                    
                    popsize = 15000000, Or_C = 3, OR_W1 = 1, OR_W2 = 2.5,
                    
                    em = 0.25, cfV0 = FALSE, cfV1 = FALSE, return_full = FALSE,
                    
                    co_inf_para = 0.00001) {
  
          
          # Génération du facteur de confusion continu C 
          
          C <- runif(n = popsize, 0.1, 3); # On commence par les noeuds racines de notre DAG
          
          # Génération du statut vaccinal V : V ~ Bernoulli(logit(10 + 0.3*C))
          
          if (cfV0 == TRUE) {
            
            V <- rep(0, popsize)
            
          } else if (cfV1 == TRUE) {
            
            V <- rep(1, popsize)
            
          } else {
            
            V <- rbinom(n = popsize, size = 1, prob = plogis(10 + 0.3*C))
            
          }
   
}


dag <- empty_dag() +
        
       node("C", type = "runif", 0.1, 3) +
       
      node("V", type = "binomial", formula = ~ 10 + 0.3*C)
  
set.seed(10)

sim_dat <- sim_from_dag(dag, n_sim = 100000)

plot(dag)

  
  