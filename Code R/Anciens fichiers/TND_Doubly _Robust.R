
# Algorithme pour l'estimation de EV en utilisant l'estimateur TNDDR

# Première étape : Diviser l'échantillon pour le cross-fitting.

# Deuxième étape : Estimer les fonctions  P_TND(V = v/ C = c, Y = 0) : score de 
                # propension parmi les témoins et P_TND(Y = 1/ V = v, C = c) en 
                # appliquant  une méthode d'apprentissage statistique. Je vais  
                # travailler  dans un premier  temps avec la méthode de la forêt 
                # d'arbres décisionnels (Random forest) pour la classification.

# Troisième étape : Estimer le RRM par : 𝜓mRR = 𝜓𝑣∕𝜓𝑣0 qui est un estimateur 
                  # à la fois doublement  robuste et efficace. (Article 03)

# Chargement des librairies nécessaires

library(sandwich)
library(randomForest)
library(dplyr)
library(ranger)
library(earth)
library(glmnet)

################################################################################

# Définir une fonction pour prédire les probabilités (Forêt aléatoire)

RandomForest <- function(dat) {
  
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
  
  ### Sur TNDdat_ctr1 
  
  # Sous-ensemble : témoins (Y == 0)
  
  TNDdat_ctr1 <- subset(TNDdat_train1, Y == 0)
  
  mod_g1_ctr <- ranger(
    
    V ~ .,
    data = subset(TNDdat_ctr1, select = -Y),
    num.trees = 200,  # Nombre d'arbres = 200
    mtry = 1,         # Set mtry to 1 for 2-3 covariates
    min.node.size = 60,  # Taille minimale d'un noeud
    sample.fraction = 0.33,
    splitrule = "extratrees",  # "Augmenter"  la part de l'aléatoire dans de la construction de l'arbre
    probability = TRUE
    
  )
  
  ### Sur TNDdat_ctr2
  
  # Sous-ensemble : témoins (Y == 0)
  
  TNDdat_ctr2 <- subset(TNDdat_train2, Y == 0)
  
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
  
  # S'assurer que newTNDdata et TNDdata ont la même structure
  
  newTNDdata_train2 <- TNDdat_train2 %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  # Prédire sur newTNDdata_train2 (ensemble autre que celui utilisé pour le premier entraînement du modèle)
  
  g1_cont[-s] <- predict(mod_g1_ctr, TNDdata = newTNDdata_train2)$predictions[, 2]
  
  newTNDdata_train1 <- TNDdat1 %>%
    mutate(V = factor(rep(1, nrow(TNDdat1)), levels = levels(TNDdat$V)))
  
  # Prédire sur newTNDdata_train2 (ensemble autre que celui utilisé pour le deuxième entraînement du modèle)
  
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
  
  mu1 <- rep(NA, nrow(TNDdat))
  mu0 <- rep(NA, nrow(TNDdat))
  
  # Prédire mu1: P(Y = 1/ V = 1) sur newTNDdata_mu1_train2
  
  newTNDdata_mu1_train2 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  mu1[-s] <- predict(Out_mu1, TNDdata = newTNDdata_mu1_train2)$predictions[, 2]
  
  # Prédire mu1: P(Y = 1/ V = 1) sur newTNDdata_mu1_train1
  
  newTNDdata_mu1_train1 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  
  mu1[s] <- predict(Out_mu2, TNDdata = newTNDdata_mu1_train1)$predictions[, 2]
  
  # Prédire mu0: P(Y = 1/ V = 0) sur newTNDdata_mu1_train2
  
  newTNDdata_mu0_train2 <- TNDdat2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  mu0[-s] <- predict(Out_mu1, TNDdata = newTNDdata_mu0_train2)$predictions[, 2]
  
  # Prédire mu0: P(Y = 1/ V = 0) sur newTNDdata_mu1_train1
  
  newTNDdata_mu0_train1 <- TNDdat1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat1)), levels = levels(TNDdat$V)))
  
  mu0[s] <- predict(Out_mu2, TNDdata = newTNDdata_mu0_train1)$predictions[, 2]
  
  # Deuxième étape : Estimer les fonctions  m0 (1 - Y ou P(Y = 0))  
  ## Entrainement du modèle de forêt aléatoire
  
  ### Sur TNDdat_ctr1   
  
  Out_m1 <- ranger(
    
    Y ~ .,
    TNDdata = subset(TNDdat1, select = -V),
    num.trees = 500,  
    mtry = 1,         
    probability = TRUE
    
  )
  
  ### Sur TNDdat_ctr2
  
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
  
  m0[-s] <- 1 - predict(Out_m1, TNDdata = select(TNDdat2, -c(V, Y)))$predictions[, 2]
  m0[s] <- 1 - predict(Out_m2, TNDdata = select(TNDdat1, -c(V, Y)))$predictions[, 2]
  
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

RandomForest <- function(dat) {
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV*dat$W2*dat$H)
 
  # Première étape : Diviser aléatoirement le jeu de données en deux parties égales 
  # Le double cross-fit
  
  set.seed(1) # pour que le résultat soit reproductible
  
  s <- sample(1:nrow(TNDdat), nrow(TNDdat) / 2)
  TNDdat_train1 <- TNDdat[s, ]
  TNDdat_train2 <- TNDdat[-s, ]
  
  # Deuxième étape : Estimer les fonctions  P_TND(V = v/ C = c, Y = 0)  
  ## Entrainement du modèle 
  
  ### Sur TNDdat_ctr1 
  
  TNDdat_ctr1 <- subset(TNDdat_train1, Y == 0)
  
  mod_g1_ctr <- cv.glmnet(
    
    V,
    X = data.matrix(subset(TNDdat_ctr1, select = -Y)),
    alpha = 1, 
    lambda = 0.2,
    standardize = TRUE
    
  )
  
  ### Sur TNDdat_ctr2
  
  # Sous-ensemble : témoins (Y == 0)
  
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
  
  # S'assurer que newTNDdata et TNDdata ont la même structure
  
  newTNDdata_train2 <- TNDdat_train2 %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  # Prédire sur newTNDdata_train2 (ensemble autre que celui utilisé pour le premier entraînement du modèle)
  
  g1_cont[-s] <- predict(mod_g1_ctr, TNDdata = newTNDdata_train2)$predictions[, 2]
  
  newTNDdata_train1 <- TNDdat1 %>%
    mutate(V = factor(rep(1, nrow(TNDdat1)), levels = levels(TNDdat$V)))
  
  # Prédire sur newTNDdata_train2 (ensemble autre que celui utilisé pour le deuxième entraînement du modèle)
  
  g1_cont[s] <- predict(mod_g2_ctr, TNDdata = newTNDdata_train1)$predictions[, 2]
  
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
  
  mu1 <- rep(NA, nrow(TNDdat))
  mu0 <- rep(NA, nrow(TNDdat))
  
  # Prédire mu1: P(Y = 1/ V = 1) sur newTNDdata_mu1_train2
  
  newTNDdata_mu1_train2 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  mu1[-s] <- predict(Out_mu1, TNDdata = newTNDdata_mu1_train2)$predictions[, 2]
  
  # Prédire mu1: P(Y = 1/ V = 1) sur newTNDdata_mu1_train1
  
  newTNDdata_mu1_train1 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  
  mu1[s] <- predict(Out_mu2, TNDdata = newTNDdata_mu1_train1)$predictions[, 2]
  
  # Prédire mu0: P(Y = 1/ V = 0) sur newTNDdata_mu1_train2
  
  newTNDdata_mu0_train2 <- TNDdat2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  mu0[-s] <- predict(Out_mu1, TNDdata = newTNDdata_mu0_train2)$predictions[, 2]
  
  # Prédire mu0: P(Y = 1/ V = 0) sur newTNDdata_mu1_train1
  
  newTNDdata_mu0_train1 <- TNDdat1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat1)), levels = levels(TNDdat$V)))
  
  mu0[s] <- predict(Out_mu2, TNDdata = newTNDdata_mu0_train1)$predictions[, 2]
  
  # Deuxième étape : Estimer les fonctions  m0 (1 - Y ou P(Y = 0))  
  ## Entrainement du modèle de forêt aléatoire
  
  ### Sur TNDdat_ctr1  
  
  Out_m1 <- cv.glmnet(
    
    Y,
    X <- data.matrix(subset(TNDdat_ctr1, select = -V)),
    alpha = 1, 
    lambda = 0.2,
    standardize = TRUE
    
  )
  
  ### Sur TNDdat_ctr2
  
  Out_m1 <- cv.glmnet(
    
    Y,
    X <- data.matrix(subset(TNDdat_ctr2, select = -V)),
    alpha = 1, 
    lambda = 0.2,
    standardize = TRUE
    
  )
  
  # Initialiser m0 pour contenir les prédictions
  
  m0 <- rep(NA, nrow(TNDdat))
  
  # Prédire mu0: P(Y = 0)
  
  m0[-s] <- 1 - predict(Out_m1, TNDdata = select(TNDdat_ctr2, -c(V, Y)))$predictions[, 2]
  m0[s] <- 1 - predict(Out_m2, TNDdata = select(TNDdat_ctr1, -c(V, Y)))$predictions[, 2]
  
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
  
  g1_cont <- dat$V
  g1_cont[-s] <- predict(mod_g1_ctr, type = "response", newdata = as.data.frame(cbind(select(TNDdat_train2, !c(V,Y)), V = rep(1, nrow(TNDdat_train2) ), Y = TNDdat_train2$Y)))
  g1_cont[s] <- predict(mod_g2_ctr, type = "response", newdata = as.data.frame(cbind(select(TNDdat_train1, !c(V,Y)), V = rep(1, nrow(TNDdat_train1)) , Y = TNDdat_train1$Y)))
  
  # Deuxième étape : Estimer les fonctions  P_TND(Y = 1/ V = v, C = c)  
  ## Entainement du modèle 
  
  ### Sur le premier ensemble d'entraînement
  
  Out_mu1 <- earth(
    
    Y ~ .,
    data = TNDdat_train1,
    glm = list(family = binomial)
    
  )
  
  ### Sur le deuxième ensemble d'entraînement
  
  Out_mu2 <- earth(
    
    Y ~ .,
    data = TNDdat_train2,
    glm = list(family = binomial)
    
  )
  
  ## Prédire les probabilités sur les ensembles tests 
  # Stockage des résultats
  
  mu1 <- TNDdat_train$Y
  mu0 <- TNDdat_train$Y
  
  mu1[-s] <- predict(Out_mu1, newdata = as.data.frame(cbind(V = 1, select(TNDdat_train2, !c(V,Y)) )), type = "response")
  mu1[s] <- predict(Out_mu2, newdata = as.data.frame(cbind(V = 1, select(TNDdat_train1, !c(V,Y)) )), type = "response")
  
  mu0[-s] <- predict(Out_mu1, newdata=as.data.frame(cbind(V = 0, select(TNDdat_train2, !c(V,Y)) )), type = "response")
  mu0[s] <- predict(Out_mu2, newdata=as.data.frame(cbind(V = 0, select(TNDdat_train1, !c(V,Y)) )), type = "response")
  
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
  
  m0 <- dat$Y
  
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

TNDDR <- function(dat){
  
}
