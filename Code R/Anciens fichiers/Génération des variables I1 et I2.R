# Charger les librairies nécessaires

install.packages("stats")
library("stats")

# On pose:

a <- 0.5/2.9
b1 <- 0.3
b2 <- 3

# Définir la fonction en b0 (en calculant l'intégrale) pour I1

f_I1 <- function (b0) {
  
  (a/b1)*(log(1 + exp(b0 + 3*b1 - b2)) + log(1 + exp(b0 + 3*b1 + b2)) - 
            
          log(1 + exp(b0 + 0.1*b1 - b2)) - log(1 + exp(b0 + 0.1*b1 + b2))) - 
    
          0.15 # Contrainte pour la probabilité marginale
  
}
uniroot(f_I1, c(-10, 10))

help("uniroot")

# Définir la fonction en b0 (en calculant l'intégrale) pour I2

# On pose:

a <- 0.5/2.9
a0 <- 0.5
a1 <- 0.3
b1 <- 0.15
b2 <- -0.1
b3 <- -3

f_I2 <- function (b0) {
  
  (a/(a1+b1))*(log(1 + exp(a0 + b0 + 3*(a1 + b1) + b2 - b3)) - log(1 + exp(a0 + b0 + 0.1*(a1 + b1) + b2 - b3)) + 
            
               log(1 + exp(a0 + b0 + 3*(a1 + b1) - b3)) - log(1 + exp(a0 + b0 + 0.1*(a1 + b1) - b3))+
              
               log(1 + exp(a0 + b0 + 3*(a1 + b1) + b2 + b3)) - log(1 + exp(a0 + b0 + 0.1*(a1 + b1) + b2 + b3)) + 
              
               log(1 + exp(a0 + b0 + 3*(a1 + b1) + b3)) - log(1 + exp(a0 + b0 + 0.1*(a1 + b1) + b3))) - 
    
               0.50 # Contrainte pour la probabilité marginale
  
}

uniroot(f_I2, c(-10, 10))

