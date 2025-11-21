# Charger les librairies nécessaires

install.packages("stats")
library("stats")
library(rje)

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

P_I2 <- function(c,b0) { # Probabilité conjointe pour I2
  
  a*(expit(b0 + b1*c + b2 - b3)*expit(a0 + a1*c) + 
                            
      expit(b0 + b1*c - b3)*expit(a0 + a1*c) + 
                            
      expit(b0 + b1*c + b2 + b3)*expit(a0 + a1*c) + 
                            
      expit(b0 + b1*c + b3)*expit(a0 + a1*c)
  )
     
}

help(integrate)
# ...	
# additional arguments to be passed to f : on utilise cet argument pour passer 
# la valeur b0 à la fonction P_I2 dans ce cas.

f_I2 <- function(b0){
  
  integrate(P_I2, 0.1, 3, b0)$value - 0.50 # Contrainte pour la probabilité marginale
  
}

uniroot(f_I2, c(-10, 10))

