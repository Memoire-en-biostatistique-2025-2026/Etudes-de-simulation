
# Algorithme pour l'estimation de EV en utilisant l'estimateur TNDDR

# Première étape : Diviser l'échantillon pour le cross-fitting.

# Deuxième étape : Estimer les fonctions  P_TND(V = v/ C = c, Y = 0) : score de 
                # propension parmi les témoins et P_TND(Y = 1/ V = v, C = c) en 
                # appliquant  une méthode d'apprentissage statistique : forêt   
                # aléatoire, régression lasso, earth_GLM et réseaux de neurones. 

# Troisième étape : Estimer le RRM par : 𝜓mRR = 𝜓𝑣∕𝜓𝑣0 qui est un estimateur 
                  # à la fois doublement  robuste et efficace. (Article 03)

# Chargement des librairies nécessaires

library(sandwich)
library(randomFoestimationst)
library(dplyr)
library(ranger)
library(earth)
library(glmnet)
library(nnet)

################################################################################

# Définir une fonction pour prédire les probabilités (Forêt aléatoire)

RandomFoest <- function(dat) {
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
  
  # Transformer les variables V et Y (variables réponses) en facteurs pour la classification
  
  TNDdat$V <- as.factor(TNDdat$V)
  TNDdat$Y <- as.factor(TNDdat$Y)
  
  # Première étape : Diviser aléatoirement le jeu de données en deux parties égales 
  # Le double cross-fit
  
  set.seed(1) # pour que le résultat soit reproductible

  s <- sample(1:nrow(TNDdat), nrow(TNDdat) / 2)
  TNDdat_train1 <- TNDdat[s, ]
  TNDdat_train2 <- TNDdat[-s, ]
  
  # Deuxième étape : Estimer les fonctions  P_TND(V = v/ C = c, Y = 0)  
  ## Entrainement du modèle de forêt aléatoire
  
  ### Premier ensemble d'entraînement : TNDdat_ctr1 
  
  TNDdat_ctr1 <- subset(TNDdat_train1, Y == 0) # Sous-ensemble : témoins (Y == 0)
  
  mod_g1_ctr <- ranger(
    
    V ~ .,
    data = subset(TNDdat_ctr1, select = -Y),
    num.trees = 200,  # Nombre d'arbestimations = 200
    mtry = 1,         # Set mtry to 1 for 2-3 covariates
    min.node.size = 60,  # Taille minimale d'un noeud
    sample.fraction = 0.33,
    splitrule = "extratrees",  # "Augmenter"  la part de l'aléatoire dans de la construction de l'arbre
    probability = TRUE
    
  )
  
  ### Deuxième ensemble d'entraînement : TNDdat_ctr2
  
  TNDdat_ctr2 <- subset(TNDdat_train2, Y == 0) # Sous-ensemble : témoins (Y == 0)
  
  mod_g2_ctr <- ranger(
    
    V ~ .,
    data = subset(TNDdat_ctr2, select = -Y),
    num.trees = 200,  
    mtry = 1,         
    min.node.size = 60,  
    sample.fraction = 0.33,
    splitrule = "extratrees",  
    probability = TRUE
    
  )
  
  ## Prédiction des probabilités sur les ensembles tests
  # Stockage des résultats
  
  g1_cont <- rep(NA, nrow(TNDdat))  
  
  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_test1 <- TNDdat_train2 %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_test2 <- TNDdat_train1 %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Prédire sur TNDdata_test1 (ensemble autre que celui utilisé pour le premier entraînement du modèle)
  
  g1_cont[-s] <- predict(mod_g1_ctr, TNDdata = newTNDdata_train2)$predictions[, 2]
  
  # Prédire sur TNDdata_test2 (ensemble autre que celui utilisé pour le deuxième entraînement du modèle)
  
  g1_cont[s] <- predict(mod_g2_ctr, TNDdata = newTNDdata_train1)$predictions[, 2]
  
  # Deuxième étape : Estimer les fonctions  P_TND(Y = 1/ V = v, C = c)  
  ## Entainement du modèle de forêt aléatoire
  
  ### Sur le premier ensemble d'entraînement
  
  Out_mu1 <- ranger(
    
    Y ~ .,
    TNDdata = TNDdat_train1,
    num.trees = 50,  
    mtry = 1,      
    probability = TRUE
    
  )
  
  ### Sur le deuxième ensemble d'entraînement
  
  Out_mu2 <- ranger(
    
    Y ~ .,
    TNDdata = TNDdat_train2,
    num.trees = 50,  
    mtry = 1,        
    probability = TRUE
    
  )
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  mu1 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 1, C = c) 

  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_mu1_test1 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_mu1_test2 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Prédire mu1: P(Y = 1/ V = 1) sur TNDdata_mu1_test1
  
  mu1[-s] <- predict(Out_mu1, TNDdata = newTNDdata_mu1_train2)$predictions[, 2]
  
  # Prédire mu1: P(Y = 1/ V = 1) sur TNDdata_mu1_test2
  
  mu1[s] <- predict(Out_mu2, TNDdata = newTNDdata_mu1_train1)$predictions[, 2]
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  mu0 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 0, C = c)
  
  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_mu0_test1 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))

  TNDdata_mu0_test2 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Prédire mu0: P(Y = 1/ V = 0) sur TNDdata_mu0_test1
  
  mu0[-s] <- predict(Out_mu1, TNDdata = TNDdata_mu0_test1)$predictions[, 2]
  
  # Prédire mu0: P(Y = 1/ V = 0) sur TNDdata_mu0_test2
  
  mu0[s] <- predict(Out_mu2, TNDdata = TNDdata_mu0_test2)$predictions[, 2]
  
  # Deuxième étape : Estimer les fonctions  m0 (1 - Y ou P(Y = 0))  
  ## Entrainement du modèle de forêt aléatoire
  
  ### Sur le premier ensemble d'entraînement   
  
  Out_m1 <- ranger(
    
    Y ~ .,
    TNDdata = subset(TNDdat1, select = -V),
    num.trees = 500,  
    mtry = 1,         
    probability = TRUE
    
  )
  
  ### Sur le deuxième ensemble d'entraînement
  
  Out_m2 <- ranger(
    
    Y ~ .,
    TNDdata = subset(TNDdat2, select = -V),
    num.trees = 500,  
    mtry = 1,         
    probability = TRUE
    
  )
  
  ## Prédiction des probabilités sur les ensembles tests
  # Stockage des résultats
  
  m0 <- rep(NA, nrow(TNDdat))
  
  m0[-s] <- 1 - predict(Out_m1, TNDdata = select(TNDdat_train2, -c(V, Y)))$predictions[, 2]
  m0[s] <- 1 - predict(Out_m2, TNDdata = select(TNDdat_train1, -c(V, Y)))$predictions[, 2]
  
  mu1 <- pmin(pmax(mu1, 0.0000001), 0.9999999)
  mu0 <- pmin(pmax(mu0, 0.0000001), 0.9999999)
  m0 <- pmin(pmax(m0, 0.0000001), 0.9999999)
  g1 <- pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  g0 <- 1 - pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1, g0 = g0, w1 = m0 / (1 - mu1), w0 = m0 / (1 - mu0)))
  
}

################################################################################
################################################################################

# Définir une fonction pour prédire les probabilités (Lasso)

Lasso <- function(dat) {
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
 
  # Première étape : Diviser aléatoirement le jeu de données en deux parties égales 
  # Le double cross-fit
  
  set.seed(1) # pour que le résultat soit reproductible
  
  s <- sample(1:nrow(TNDdat), nrow(TNDdat) / 2)
  TNDdat_train1 <- TNDdat[s, ]
  TNDdat_train2 <- TNDdat[-s, ]
  
  # Deuxième étape : Estimer les fonctions  P_TND(V = v/ C = c, Y = 0)  
  ## Entrainement du modèle 
  
  ### Premier ensemble d'entraînement : TNDdat_ctr1  
  
  TNDdat_ctr1 <- subset(TNDdat_train1, Y == 0)
  
  mod_g1_ctr <- cv.glmnet(
    
    V,
    X = data.matrix(subset(TNDdat_ctr1, select = -Y)),
    alpha = 1, 
    lambda = 0.2,
    standardize = TRUE
    
  )
  
  ### Deuxième ensemble d'entraînement : TNDdat_ctr2
  
  TNDdat_ctr2 <- subset(TNDdat_train2, Y == 0)
  
  mod_g1_ctr <- cv.glmnet(
    
    V,
    X = data.matrix(subset(TNDdat_ctr2, select = -Y)),
    alpha = 1, 
    lambda = 0.2,
    standardize = TRUE
    
  )
  
  ## Prédiction des probabilités sur les ensembles tests
  # Stockage des résultats
  
  g1_cont <- rep(NA, nrow(TNDdat))  
  
  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_test1 <- data.matrix(as.data.frame(cbind(select(TNDdat_train2, !c(V,Y)),
                                                    V = rep(1, nrow(TNDdat_train2) ), 
                                                    Y = TNDdat_train2$Y)))
  
  TNDdata_test2 <- data.matrix(as.data.frame(cbind(select(TNDdat_train2, !c(V,Y)),
                                                   V = rep(1, nrow(TNDdat_train2) ), 
                                                   Y = TNDdat_train2$Y)))
  
  # Prédire sur TNDdata_test1 (ensemble autre que celui utilisé pour le premier entraînement du modèle)
  
  g1_cont[-s] <- predict(mod_g1_ctr, newx = TNDdata_test1)
  
  # Prédire sur TNDdata_test2 (ensemble autre que celui utilisé pour le deuxième entraînement du modèle)
  
  g1_cont[s] <- predict(mod_g2_ctr, newx = TNDdata_test2)
  
  # Deuxième étape : Estimer les fonctions  P_TND(Y = 1/ V = v, C = c)  
  ## Entainement du modèle de forêt aléatoire
  
  ### Sur TNDdat_train1
  
  Out_mu1 <- cv.glmnet(
    
    Y,
    X = data.matrix(TNDdat_train1),
    alpha = 1, 
    lambda = 0.2,
    standardize = TRUE
    
  )
  
  ### Sur TNDdat_train2
  
  Out_mu1 <- cv.glmnet(
    
    Y,
    X = data.matrix(TNDdat_train1),
    alpha = 1, 
    lambda = 0.2,
    standardize = TRUE
    
  )
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  mu1 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 1, C = c) 
  
  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_mu1_test1 <- data.matrix(as.data.frame(cbind(V = 1, select(TNDdat_train2, !c(V,Y)))))
  
  TNDdata_mu1_test2 <- data.matrix(as.data.frame(cbind(V = 1, select(TNDdat_train1, !c(V,Y)))))
  
  # Prédire mu1: P(Y = 1/ V = 1) sur TNDdata_mu1_test1
  
  mu1[s] <- predict(Out_mu2, newx = TNDdata_mu1_test1)
  
  # Prédire mu1: P(Y = 1/ V = 1) sur newTNDdata_mu1_train2
  
  mu1[-s] <- predict(Out_mu1, newx = TNDdata_mu1_test2)
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  mu0 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 0, C = c) 
  
  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_mu1_test1 <- data.matrix(as.data.frame(cbind(V = 0, select(TNDdat_train2, !c(V,Y)))))
  
  TNDdata_mu1_test2 <- data.matrix(as.data.frame(cbind(V = 0, select(TNDdat_train1, !c(V,Y)))))
  
  # Prédire mu0: P(Y = 1/ V = 0) sur TNDdata_mu1_test1
  
  mu0[s] <- predict(Out_mu2, newx = TNDdata_mu1_test1)
  
  # Prédire mu1: P(Y = 1/ V = 0) sur newTNDdata_mu1_train2
  
  mu0[-s] <- predict(Out_mu1, newx = TNDdata_mu1_test2)
  
  # Deuxième étape : Estimer les fonctions  m0 (1 - Y ou P(Y = 0))   
  ## Entainement du modèle 
  
  ### Sur le premier ensemble d'entraînement  
  
  Out_m1 <- cv.glmnet(
    
    Y,
    X <- data.matrix(subset(TNDdat_ctr1, select = -V)),
    alpha = 1, 
    lambda = 0.2,
    standardize = TRUE
    
  )
  
  ### Sur le deuxième ensemble d'entraînement
  
  Out_m1 <- cv.glmnet(
    
    Y,
    X <- data.matrix(subset(TNDdat_ctr2, select = -V)),
    alpha = 1, 
    lambda = 0.2,
    standardize = TRUE
    
  )
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  m0 <- rep(NA, nrow(TNDdat))
  
  m0[-s] <- 1 - predict(Out_m1, newx = data.matrix(select(TNDdat_train2, !c(V,Y)), type = "response"))
  m0[s] <- 1 - predict(Out_m2, newx = data.matrix(select(TNDdat_train1, !c(V,Y)), type = "response"))
  
  mu1 <- pmin(pmax(mu1, 0.0000001), 0.9999999)
  mu0 <- pmin(pmax(mu0, 0.0000001), 0.9999999)
  m0 <- pmin(pmax(m0, 0.0000001), 0.9999999)
  g1 <- pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  g0 <- 1 - pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1, g0 = g0, w1 = m0 / (1 - mu1), w0 = m0 / (1 - mu0)))
  
}



################################################################################
################################################################################

# Définir une fonction pour prédire les probabilités (earth_GLM)

Mars <- function(dat){
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
  
  # Première étape : Diviser aléatoirement le jeu de données en deux parties égales 
  # Le double cross-fit
  
  set.seed(1) # pour que le résultat soit reproductible
  
  s <- sample(1:nrow(TNDdat), nrow(TNDdat) / 2)
  
  TNDdat_train1 <- TNDdat[s, ]
  TNDdat_train2 <- TNDdat[-s, ]
  
  # Deuxième étape : Estimer les fonctions  P_TND(V = v/ C = c, Y = 0)  
  ## Entrainement du modèle de forêt aléatoire
  
  ### Sur le premier ensemble d'entraînement
  
  TNDdat_train_ctr1 <- subset(TNDdat_train1, Y == 0)
  
  mod_g1_ctr <- earth(
    
    V ~ .,
    data = subset(TNDdat_train_ctr1, select = -Y),
    glm = list(family = binomial)
    
  )
  
  ### Sur le deuxième ensemble d'entraînement
  
  TNDdat_train_ctr2 <- subset(TNDdat_train2, Y==0)
  
  mod_g2_ctr <- earth(
    
    V ~ .,
    data = subset(TNDdat_train_ctr2, select = -Y),
    glm = list(family = binomial)
    
  )
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  g1_cont <- rep(NA, nrow(TNDdat)) 
  
  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_test1 <-  newdata = as.data.frame(cbind(select(TNDdat_train2, !c(V,Y)),
                                                  V = rep(1, nrow(TNDdat_train2) ), 
                                                  Y = TNDdat_train2$Y))
  
  TNDdata_test2 <- newdata = as.data.frame(cbind(select(TNDdat_train1, !c(V,Y)),
                                                 V = rep(1, nrow(TNDdat_train1) ), 
                                                 Y = TNDdat_train1$Y))
  
  # Prédire sur TNDdata_test1 (ensemble autre que celui utilisé pour le premier entraînement du modèle)
  
  g1_cont[-s] <- predict(mod_g1_ctr, type = "response", newdata = TNDdata_test1)
  
  # Prédire sur TNDdata_test2 (ensemble autre que celui utilisé pour le deuxième entraînement du modèle)
  
  g1_cont[s] <- predict(mod_g2_ctr, type = "response", newdata = TNDdata_test2)
  
  # Deuxième étape : Estimer les fonctions  P_TND(Y = 1/ V = v, C = c)  
  ## Entainement du modèle 
  
  ### Sur le premier ensemble d'entraînement
  
  Out_mu1 <- earth(
    
    Y ~ .,
    data = TNDdat_train1,
    glm = list(family = binomial),
    degree = 2
    
  )
  
  ### Sur le deuxième ensemble d'entraînement
  
  Out_mu2 <- earth(
    
    Y ~ .,
    data = TNDdat_train2,
    glm = list(family = binomial),
    degree = 2
    
  )
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  mu1 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 1, C = c)
  
  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_mu1_test1 <- as.data.frame(cbind(V = 1, select(TNDdat_train2, !c(V,Y)) ))
  
  TNDdata_mu1_test2 <- as.data.frame(cbind(V = 1, select(TNDdat_train1, !c(V,Y)) ))
  
  # Prédire mu1: P(Y = 1/ V = 1) sur TNDdata_mu1_test1
  
  mu1[-s] <- predict(Out_mu1, newdata = TNDdata_mu1_test1, type = "response")
  
  # Prédire mu1: P(Y = 1/ V = 1) sur TNDdata_mu1_test2
  
  mu1[s] <- predict(Out_mu2, newdata = TNDdata_mu1_test2, type = "response")
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  mu0 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 0, C = c)
  
  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_mu1_test1 <- as.data.frame(cbind(V = 0, select(TNDdat_train2, !c(V,Y)) ))
  
  TNDdata_mu1_test2 <- as.data.frame(cbind(V = 0, select(TNDdat_train1, !c(V,Y)) ))
  
  # Prédire mu0: P(Y = 1/ V = 0) sur TNDdata_mu1_test1
  
  mu0[-s] <- predict(Out_mu1, newdata = TNDdata_mu1_test1, type = "response")
  
  # Prédire mu0: P(Y = 1/ V = 0) sur TNDdata_mu1_test2
  
  mu0[s] <- predict(Out_mu2, newdata = TNDdata_mu1_test2, type = "response")
  
  # Deuxième étape : Estimer les fonctions  m0 (1 - Y ou P(Y = 0))   
  ## Entainement du modèle 
  
  ### Sur le premier ensemble d'entraînement

  Out_m1 <- earth(
    
    Y ~ .,
    data = subset(TNDdat_train1, select = -V),
    glm = list(family = binomial)
    
  )
  
  ### Sur le deuxième ensemble d'entraînement
  
  Out_m2 <- earth(
    
    Y ~ .,
    data = subset(TNDdat_train2, select = -V),
    glm = list(family = binomial)
    
  )
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  m0 <- rep(NA, nrow(TNDdat))
  
  m0[-s] <- 1 - predict(Out_m1, newdata = select(TNDdat_train2, !c(V,Y)), type = "response")
  m0[s] <- 1 - predict(Out_m2, newdata = select(TNDdat_train1, !c(V,Y)), type = "response")
  
  mu1 <- pmin(pmax(mu1, 0.0000001), 0.9999999)
  mu0 <- pmin(pmax(mu0, 0.0000001), 0.9999999)
  m0 <- pmin(pmax(m0, 0.0000001), 0.9999999)
  g1 <- pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  g0 <- 1 - pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1,g0 = g0, w1 = m0 / (1 - mu1),w0 = m0 / (1 - mu0)))
  
}

################################################################################
################################################################################

RN <- function(dat) {
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
  
  # Transformer les variables V et Y (variables réponses) en facteurs pour la classification
  
  TNDdat$V <- as.factor(TNDdat$V)
  TNDdat$Y <- as.factor(TNDdat$Y)
  
  # Première étape : Diviser aléatoirement le jeu de données en deux parties égales 
  # Le double cross-fit
  
  set.seed(1) # pour que le résultat soit reproductible
  
  s <- sample(1:nrow(TNDdat), nrow(TNDdat) / 2)
  TNDdat_train1 <- TNDdat[s, ]
  TNDdat_train2 <- TNDdat[-s, ]
  
  # Deuxième étape : Estimer les fonctions  P_TND(V = v/ C = c, Y = 0)  
  ## Entrainement du modèle de forêt aléatoire
  
  ### Premier ensemble d'entraînement : TNDdat_ctr1 
  
  TNDdat_ctr1 <- subset(TNDdat_train1, Y == 0) # Sous-ensemble : témoins (Y == 0)
  
  mod_g1_ctr <- nnet( # La fonction d’activation par défaut est la fonction sigmoïde
                      # Cela aide le réseau à introduire la non-linéarité
    V ~ .,
    data = subset(TNDdat_ctr1, select = -Y),
    size = 5, # Le nombre de nœuds dans la couche cachée
    maxit = 50 # Le paramètre fixe le nombre maximal d’itérations pour l'e processus d’entraînement'entraînement
    
  )
  
  ### Deuxième ensemble d'entraînement : TNDdat_ctr2
  
  TNDdat_ctr2 <- subset(TNDdat_train2, Y == 0) # Sous-ensemble : témoins (Y == 0)
  
  mod_g2_ctr <- nnet(
    
    V ~ .,
    data = subset(TNDdat_ctr2, select = -Y),
    size = 5, 
    maxit = 50
    
  )
  
  ## Prédiction des probabilités sur les ensembles tests
  # Stockage des résultats
  
  g1_cont <- rep(NA, nrow(TNDdat))  
  
  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_test1 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_test2 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Prédire sur TNDdata_test1 (ensemble autre que celui utilisé pour le premier entraînement du modèle)
  
  g1_cont[-s] <- predict(mod_g1_ctr, TNDdata = newTNDdata_train2, type = "raw")
  
  # Prédire sur TNDdata_test2 (ensemble autre que celui utilisé pour le deuxième entraînement du modèle)
  
  g1_cont[s] <- predict(mod_g2_ctr, TNDdata = newTNDdata_train1, type = "raw")
  
  # Deuxième étape : Estimer les fonctions  P_TND(Y = 1/ V = v, C = c)  
  ## Entainement du modèle de forêt aléatoire
  
  ### Sur le premier ensemble d'entraînement
  
  Out_mu1 <- nnet(
    
    Y ~ .,
    data = TNDdat_train1,
    size = 5, 
    maxit = 50
    
  )
  
  ### Sur le deuxième ensemble d'entraînement
  
  Out_mu2 <- nnet(
    
    Y ~ .,
    data = TNDdat_train2,
    size = 5, maxit = 50
    
  )
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  mu1 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 1, C = c) 
  
  # S'assurer à chaque fois que les ensemble d'entraînemnt et de test ont la même structure
  
  TNDdata_mu1_test1 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_mu1_test2 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Prédire mu1: P(Y = 1/ V = 1) sur TNDdata_mu1_test1
  
  mu1[-s] <- predict(Out_mu1, TNDdata = newTNDdata_mu1_train2, type = "raw")
  
  # Prédire mu1: P(Y = 1/ V = 1) sur TNDdata_mu1_test2
  
  mu1[s] <- predict(Out_mu2, TNDdata = newTNDdata_mu1_train1, type = "raw")
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  mu0 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 0, C = c)
  
  # S'assurer à chaque fois que les ensembles d'entraînemnt et de test ont la même structure
  
  TNDdata_mu0_test1 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_mu0_test2 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Prédire mu0: P(Y = 1/ V = 0) sur TNDdata_mu0_test1
  
  mu0[-s] <- predict(Out_mu1, TNDdata = TNDdata_mu0_test1, type = "raw")
  
  # Prédire mu0: P(Y = 1/ V = 0) sur TNDdata_mu0_test2
  
  mu0[s] <- predict(Out_mu2, TNDdata = TNDdata_mu0_test2, type = "raw")
  
  # Deuxième étape : Estimer les fonctions  m0 (1 - Y ou P(Y = 0))  
  ## Entrainement du modèle de forêt aléatoire
  
  ### Sur le premier ensemble d'entraînement   
  
  Out_m1 <- nnet(
    
    Y ~ .,
    data = subset(TNDdat_train1, select = -V),
    size = 5,
    maxit = 50
    
  )
  
  ### Sur le deuxième ensemble d'entraînement
  
  Out_m2 <- nnet(
    
    Y ~ .,
    data = subset(TNDdat_train2, select = -V),
    size = 5, 
    maxit = 50
    
  )
  
  ## Prédiction des probabilités sur les ensembles tests
  # Stockage des résultats
  
  m0 <- rep(NA, nrow(TNDdat))
  
  m0[-s] <- 1 - predict(Out_m1, TNDdata = select(TNDdat_train2, -c(V, Y)), type = "raw")
  m0[s] <- 1 - predict(Out_m2, TNDdata = select(TNDdat_train1, -c(V, Y)), type = "raw")
  
  mu1 <- pmin(pmax(mu1, 0.0000001), 0.9999999)
  mu0 <- pmin(pmax(mu0, 0.0000001), 0.9999999)
  m0 <- pmin(pmax(m0, 0.0000001), 0.9999999)
  g1 <- pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  g0 <- 1 - pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1, g0 = g0, w1 = m0 / (1 - mu1), w0 = m0 / (1 - mu0)))
  
  
}
  
################################################################################
################################################################################

# Troisième étape : Estimer le RRM par : 𝜓mRR = 𝜓𝑣∕𝜓𝑣0 qui est un estimateur 
# à la fois doublement  robuste et efficace. (Article 03)
# Formules mathématiques utilisées : voir Article 03 : A Double Machine Learning Approach...

# EIF (the efficient influence function) : La fonction d’influence empirique est une mesure 
# de la dépendance de l’estimateur à la valeur de l’un des points de l’échantillon. C’est 
# une mesure sans modèle dans le sens où elle repose simplement sur le calcul e l’estimateur 
# à nouveau avec un autre échantillon.

# D'abord avec la forêt aléatoire

TNDDR <- function(dat){
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
  estimations <- RandomFoest(dat)
  
  # Estimation du 𝜓𝑣: 𝜓𝑣 = 𝜓𝑣(ℙTND)=𝔼TND[𝜇𝑣(𝒄)*𝜔𝑣(𝒄)] avec
  #  𝜇𝑣(𝒄) = ℙTND (𝑌 = 1|𝑉 = 1,𝑪 = 𝒄),
  #  𝜔𝑣(𝒄) = 𝜋𝑣(𝒄) / 𝜋0𝑣(𝒄)
  
  A.1 <- ((1 - TNDdat$Y)*(TNDdat$V - estimations$g1))/(estimations$g1* (1 -estimations$mu1))
  psi.1 <- mean(TNDdat$Y*TNDdat$V/estimations$g1 - estimations$mu1*A.1)
  
  # Estimation du 𝜓𝑣0 : 𝜓𝑣0(ℙTND)=𝔼TND[𝜇𝑣0(𝒄)*𝜔𝑣0(𝒄)] avec
  #  𝜇𝑣0(𝒄) = ℙTND (𝑌 = 1|𝑉 = 0,𝑪 = 𝒄),

  A.0 <- ((1 - TNDdat$Y)*((1-TNDdat$V) - estimations$g0))/(estimations$g0* (1 - estimations$mu0))
  psi.0 <- mean(TNDdat$Y*(1-TNDdat$V)/estimations$g0 - estimations$mu0*A.0)
  
  # Estimation du RRM par : 𝜓mRR = 𝜓𝑣∕𝜓𝑣0
  
  RRm <- pmin(pmax(psi.1/psi.0, 0.001), 0.999)
  
  # Intervalles de confiance
  
  ## Méthode 01: Pour l'intervalle de confiance de TNDDR, on peut utiliser une transformation 
  ## logarithmique pour améliorer la précision de l'approximation normale.
  
  ## Estimation de var(ln(𝜓eif mRR)) = var(𝔼𝕀𝔽(ln(𝜓𝑣∕𝜓𝑣0 )))
  
  log_eif_RRm <-  ((TNDdat$Y*TNDdat$V/estimations$g1 - estimations$mu1*A.1 - psi.1)/psi.1) - ((TNDdat$Y*(1-TNDdat$V)/estimations$g0 - estimations$mu0*A.0 - psi.0)/ psi.0)
  var_log_eif_RRm <-  var(log_eif_RRm)/nrow(TNDdat)
  
  ## Premier intervalle de confiance pour le RRm
  
  IC_inf1 <- exp(log(RRm) - 1.96 * sqrt(var_log_eif_RRm) )
  IC_sup1 <- exp(log(RRm) + 1.96 * sqrt(var_log_eif_RRm) )
  
  eifpsi <- (TNDdat$Y*TNDdat$V/estimations$g1 - estimations$mu1*A.1 - psi.1)/psi.0 - RRm*(TNDdat$Y*(1-TNDdat$V)/estimations$g0 - estimations$mu0*A.0 - psi.0)/psi.0
  var <- var(eifpsi)/nrow(TNDdat)
  
  ## Deuxième intervalle de confiance pour le RRm
  
  IC_inf2 <- RRm - 1.96 * sqrt(var)
  IC_sup2 <- RRm + 1.96 * sqrt(var)
  
  varn2 <- mean((RRm* (TNDdat$Y*(1-TNDdat$V)/estimations$g0 - estimations$mu0*A.0) - (TNDdat$Y*TNDdat$V/estimations$g1 - estimations$mu1*A.1) )^2)
  denJ <- psi.0^2
  var2 <- varn2/(denJ * nrow(TNDdat))
  
  ## Troisième intervalle de confiance pour le RRm
  
  IC_inf3 <- RRm - 1.96 * sqrt(var2)
  IC_sup3 <- RRm + 1.96 * sqrt(var2)
  
  return(list( est = RRm, varln = var_log_eif_RRm, var = var, IC1 =  c(IC_inf1,IC_sup1), IC2 =  c(IC_inf2,IC_sup2), CI3 =  c(IC_inf3,IC_sup3)))
  
}
