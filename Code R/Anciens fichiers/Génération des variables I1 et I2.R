# Charger les librairies nécessaires

install.packages("stats")
library("stats")

# On pose:

a <- 0.5/2.9
b1 <- 0.3
b2 <- 0.1

# Définir la fonction en b0 (en calculant l'intégrale)

f <- function (b0) {
  
  (a/b1)*(log(1 + exp(b0 + 3*b1 - b2)) + log(1 + exp(b0 + 3*b1 + b2)) - 
            
          log(1 + exp(b0 + 0.1*b1 - b2)) - log(1 + exp(b0 + 0.1*b1 + b2))) - 
    
          0.15 # Contrainte pour la probabilité marginale
  
}
uniroot(f, c(-10, 10))

help("uniroot")
