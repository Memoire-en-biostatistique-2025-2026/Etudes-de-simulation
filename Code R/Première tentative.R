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
  
  # Génération du statut vaccinal V : V ~ Bernoulli(logit(a0 + a1*C))
  
  V <- rbinom(n = popsize, size = 1, prob = plogis(0.5 + 0.3*C))
    
  # Génération des infections I1 : I1 ~ Bernoulli(logit(a0 + a1*C)) et 
  #                           I2 : I2 ~ Bernoulli(logit(a0 + a1*C + a2*V))
  
  
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(-2.5 + 0.35*C))
  
  I2 <- rbinom(n = popsize, size = 1, prob = plogis(0 + 0.15*C - 0.1*V))
  
  # Calcul du pourcentage des co-infections
  
  co_inf <- sum(I1 == 1 & I2 == 1)
  
  per_co_inf <- co_inf / popsize * 100
  # 7.120667
  
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
    
  dat <- as.data.frame(cbind(Infec_RSV = I2, Infec = I1, H = H, W = W, V = V,
                               C = C)) # Virus respiratoire syncytial RSV
  
  return(dat)
  
}  

dat <- datagen(ssize = 5000)

################################################################################ 

datagen.cont <- function(seed = sample(1:1000000, size = 1), popsize = 150000) {
  
  
  # Génération du facteur de confusion continu C 
  
  C <- runif(n = popsize, 0.1, 3); # On commence par les noeuds racines de notre DAG
  
  # Génération des infections I1 : I1 ~ Bernoulli(logit(a0 + a1*C)) et 
  #                           I2 : I2 ~ Bernoulli(logit(a0 + a1*C + a2*V))
  
  
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(-2.5 + 0.35*C))
  
  I2_0 <- rbinom(n = popsize, size = 1, prob = plogis(0 + 0.15*C - 0.1*0))
  I2_1 <- rbinom(n = popsize, size = 1, prob = plogis(0 + 0.15*C - 0.1*1))
  
  # Calcul du pourcentage des co-infections
  
  co_inf_0 <- sum(I1 == 1 & I2_0 == 1)
  
  per_co_inf_0 <- co_inf_0 / popsize * 100
  # 7.120667
  
  co_inf_1 <- sum(I1 == 1 & I2_1 == 1)
  
  per_co_inf_1 <- co_inf_1 / popsize * 100
  # 7.120667
  
  
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
  
  return(dat)
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
################################################################################

#### Calcul des vraies valeurs ####

## Regression logistique ##

dat <- datagen.cont(seed = 94178, popsize = 10000000)

dat0 <- data.frame(C = dat$C, V = 0, Y = dat$I2_0*dat$W2_0*dat$H_0)

dat1 <- data.frame(C = dat$C, V = 1, Y = dat$I2_1*dat$W2_1*dat$H_1) 

dat_complet <- rbind(dat0, dat1)

vraiRRc <- glm(Y ~ V + C, family = binomial(link = "log"), data = dat_complet)

vrai.EV.logistique <- 1 - exp(coef(vraiRRc)[2])
# 0.3624796

## Autres methodes ##

vraiRRm <- mean(dat1$Y)/mean(dat0$Y)

vrai.EV.autres <- 1 - vraiRRm
# 0.3665348
################################################################################

#### Paramètres de la simulation ####

## Regression logistique ##

set.seed(1) # Pour avoir toujours les mêmes germes
nrep <- 10

seeds_list <- sample(1:1000000, size = nrep)

l_vraiRRc <- rep(NA, nrep)
vrai.EV.autres <- rep(NA, nrep)

j = 1;
for (i in seeds_list) {
  
  # for (j in 1:50){

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
  
vrai.EV.autres
# [1] 0.3639892 0.3659440 0.3674136 0.3673475 0.3680050 0.3652723 0.3669136 0.3654393 0.3664639
# [10] 0.3673925
mean(vrai.EV.autres)
# 0.3664181
sd(vrai.EV.autres)
# 0.001245681

