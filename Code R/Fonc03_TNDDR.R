
# Install and load required packages

library(sandwich)
library(randomForest)
library(dplyr)
library(ranger)
library(earth)
library(glmnet)
library(nnet)

# Define different functions to predict probabilities

## Random forest

RandomForest <- function(dat) {
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV)
  
  # Convert the variables V and Y (response variables) into factors for classification
  
  TNDdat$V <- as.factor(TNDdat$V)
  TNDdat$Y <- as.factor(TNDdat$Y)
  
  # Step 1: Randomly split the dataset into two equal parts 
  # Double CrossFit
  
  s <- sample(1:nrow(TNDdat), nrow(TNDdat) / 2)
  TNDdat_train1 <- TNDdat[s, ]
  TNDdat_train2 <- TNDdat[-s, ]
  
  # Step 2: Estimate the probability  P_TND(V = v/ C = c, Y = 0)  
  ## Training the random forest model
  
  ### First training set: TNDdat_ctr1 
  
  TNDdat_ctr1 <- subset(TNDdat_train1, Y == 0) # Subset: controls (Y == 0)
  
  mod_g1_ctr <- ranger(
    
    V ~ .,
    data = subset(TNDdat_ctr1, select = -Y),
    num.trees = 500,  # Number of trees in the random forest
    mtry = 1,         # Candidate partitioning variables for each partition: Set mtry to 1 for 2–3 covariates
    min.node.size = 49,  # Minimum node size
    sample.fraction = 0.8,
    splitrule = "extratrees",  # “Increase” the degree of randomness in tree construction
    probability = TRUE
    
  )
  
  ### Second training set: TNDdat_ctr2
  
  TNDdat_ctr2 <- subset(TNDdat_train2, Y == 0) # Subset: controls (Y == 0)
  
  mod_g2_ctr <- ranger(
    
    V ~ .,
    data = subset(TNDdat_ctr2, select = -Y),
    num.trees = 500,  
    mtry = 1,         
    min.node.size = 49,  
    sample.fraction = 0.8,
    splitrule = "extratrees",  
    probability = TRUE
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  g1_cont <- rep(NA, nrow(TNDdat))  
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_test1 <- TNDdat_train2 %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_test2 <- TNDdat_train1 %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Predict on TNDdata_test1 (a dataset other than the one used for the model's initial training)
  
  g1_cont[-s] <- predict(mod_g1_ctr, data = TNDdata_test1)$predictions[, 2]
  
  # Predict on TNDdata_test2 (a dataset other than the one used for the second model training)
  
  g1_cont[s] <- predict(mod_g2_ctr, data = TNDdata_test2)$predictions[, 2]
  
  # Step 3: Estimate the probability  P_TND(Y = 1/ V = v, C = c)  
  ## Training the random forest model
  
  ### First training set
  
  Out_mu1 <- ranger(
    
    Y ~ .,
    data = TNDdat_train1,
    num.trees = 50,  
    mtry = 1,
    min.node.size = 11,
    probability = TRUE
    
  )
  
  ### Second training set
  
  Out_mu2 <- ranger(
    
    Y ~ .,
    data = TNDdat_train2,
    num.trees = 50,  
    mtry = 1,   
    min.node.size = 11,
    probability = TRUE
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  mu1 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 1, C = c) 
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_mu1_test1 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_mu1_test2 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Predict mu1: P(Y = 1/ V = 1) on TNDdata_mu1_test1
  
  mu1[-s] <- predict(Out_mu1, data = TNDdata_mu1_test1)$predictions[, 2]
  
  # Predict mu1: P(Y = 1/ V = 1) on TNDdata_mu1_test2
  
  mu1[s] <- predict(Out_mu2, data = TNDdata_mu1_test2)$predictions[, 2]
  
  ## Predicting probabilities on test sets
  # Store results
  
  mu0 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 0, C = c)
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_mu0_test1 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_mu0_test2 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Predict mu0: P(Y = 1/ V = 0) on TNDdata_mu0_test1
  
  mu0[-s] <- predict(Out_mu1, data = TNDdata_mu0_test1)$predictions[, 2]
  
  # Predict mu0: P(Y = 1/ V = 0) on TNDdata_mu0_test2
  
  mu0[s] <- predict(Out_mu2, data = TNDdata_mu0_test2)$predictions[, 2]
  
  # Step 4: Estimate the probability  m0 (1 - Y or P(Y = 0))  
  ## Training the random forest model
  
  ### First training set  
  
  Out_m1 <- ranger(
    
    Y ~ .,
    data = subset(TNDdat_train1, select = -V),
    num.trees = 500,  
    mtry = 1,  
    min.node.size = 49,
    probability = TRUE
    
  )
  
  ### Second training set
  
  Out_m2 <- ranger(
    
    Y ~ .,
    data = subset(TNDdat_train2, select = -V),
    num.trees = 500,  
    mtry = 1,  
    min.node.size = 49,
    probability = TRUE
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  m0 <- rep(NA, nrow(TNDdat))
  
  m0[-s] <- 1 - predict(Out_m1, data = select(TNDdata_test1, -c(V, Y)))$predictions[, 2]
  m0[s] <- 1 - predict(Out_m2, data = select(TNDdata_test2, -c(V, Y)))$predictions[, 2]
  
  mu1 <- pmin(pmax(mu1, 0.01), 0.99)
  mu0 <- pmin(pmax(mu0, 0.01), 0.99)
  m0 <- pmin(pmax(m0, 0.01), 0.99)
  g1 <- pmin(pmax(g1_cont, 0.01), 0.99)
  g0 <- 1 - pmin(pmax(g1_cont, 0.01), 0.99)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1, g0 = g0, w1 = m0 / (1 - mu1), w0 = m0 / (1 - mu0)))
  
}

################################################################################

# GLM

PM <- function(dat) {
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV)
  
  # Step 1: Randomly split the dataset into two equal parts 
  # Double CrossFit
  
  s <- sample(1:nrow(TNDdat), nrow(TNDdat) / 2)
  TNDdat_train1 <- TNDdat[s, ]
  TNDdat_train2 <- TNDdat[-s, ]
  
  # Step 2: Estimate the probability P_TND(V = v/ C = c, Y = 0)  
  ## Training the Lasso model
  
  ### First training set: TNDdat_ctr1 
  
  TNDdat_ctr1 <- subset(TNDdat_train1, Y == 0)
  
  mod_g1_ctr <- glm( 
    
    V ~ C,
    data = TNDdat_ctr1,
    family = binomial(link = "logit")
    
  )
  
  ### Second training set : TNDdat_ctr2
  
  TNDdat_ctr2 <- subset(TNDdat_train2, Y == 0)
  
  mod_g2_ctr <- glm( 
    
    V ~ C,
    data = TNDdat_ctr2,
    family = binomial(link = "logit")
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  g1_cont <- rep(NA, nrow(TNDdat))  
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_test1 <- data.matrix(as.data.frame(cbind(select(TNDdat_train2, !c(V))) 
  ))
  
  
  TNDdata_test2 <- data.matrix(as.data.frame(cbind(select(TNDdat_train1, !c(V)) 
  )))
  
  # Predict on TNDdata_test1 (a dataset other than the one used for the model's initial training)
  
  g1_cont[-s] <- predict(mod_g1_ctr, newdata = as.data.frame(TNDdata_test1), type = "response")
  
  # Predict on TNDdata_test2 (a dataset other than the one used for the second model training)
  
  g1_cont[s] <- predict(mod_g2_ctr, newdata = as.data.frame(TNDdata_test2), type = "response")
  
  # Step 3: Estimate the probability  P_TND(Y = 1/ V = v, C = c)  
  ## Training the random forest model
  
  ### First training set
  
  Out_mu1 <- glm( 
    
    Y ~ C + V,
    data = TNDdat_train1,
    family = binomial(link = "logit")
    
  )
    
  ### Second training set
  
  Out_mu2 <- glm( 
      
      Y ~ C + V,
      data = TNDdat_train2,
      family = binomial(link = "logit")
      
    )
    
  ## Predicting probabilities on test sets
  # Store results
  
  mu1 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 1, C = c) 
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_mu1_test1 <- data.matrix(as.data.frame(cbind(V = 1, select(TNDdat_train2, !c(Y, V)))))
  
  TNDdata_mu1_test2 <- data.matrix(as.data.frame(cbind(V = 1, select(TNDdat_train1, !c(Y, V)))))
  
  # Predict mu1: P(Y = 1/ V = 1) on TNDdata_mu1_test1
  
  mu1[-s] <- predict(Out_mu1, newdata = as.data.frame(TNDdata_mu1_test1), type = "response")
  
  # Predict mu1: P(Y = 1/ V = 1) on TNDdata_mu1_test2
  
  mu1[s] <- predict(Out_mu2, newdata = as.data.frame(TNDdata_mu1_test2), type = "response")
  
  ## Predicting probabilities on test sets
  # Store results
  
  mu0 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 0, C = c) 
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_mu1_test1 <- data.matrix(as.data.frame(cbind(V = 0, select(TNDdat_train2, !c(Y, V)))))
  
  TNDdata_mu1_test2 <- data.matrix(as.data.frame(cbind(V = 0, select(TNDdat_train1, !c(Y, V)))))
  
  # Predict mu0: P(Y = 1/ V = 0) on TNDdata_mu1_test1
  
  mu0[-s] <- predict(Out_mu1, newdata = as.data.frame(TNDdata_mu1_test1), type = "response")
  
  # Predict mu0: P(Y = 1/ V = 0) on TNDdata_mu1_test2
  
  mu0[s] <- predict(Out_mu2, newdata = as.data.frame(TNDdata_mu1_test2), type = "response")
  
  # Step 4: Estimate the probability  m0 (1 - Y or P(Y = 0)) 
  ## Training the Lasso model
  
  ### First training set  
  
  Out_m1 <- glm( 
    
    Y ~ C,
    data = TNDdat_train1,
    family = binomial(link = "logit")
    
  )
  
  ### Second training set
  
  Out_m2 <- glm( 
    
    Y ~ C,
    data = TNDdat_train2,
    family = binomial(link = "logit")
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  m0 <- rep(NA, nrow(TNDdat))
  
  m0[-s] <- 1 - predict(Out_m1, newdata = as.data.frame(cbind(data.matrix(select(TNDdat_train2, !c(V,Y))), 0)), type = "response")
  m0[s] <- 1 - predict(Out_m2, newdata = as.data.frame(cbind(data.matrix(select(TNDdat_train1, !c(V,Y))), 0)), type = "response")
  
  mu1 <- pmin(pmax(mu1, 0.01), 0.99)
  mu0 <- pmin(pmax(mu0, 0.01), 0.99)
  m0 <- pmin(pmax(m0, 0.01), 0.99)
  g1 <- pmin(pmax(g1_cont, 0.01), 0.99)
  g0 <- 1 - pmin(pmax(g1_cont, 0.01), 0.99)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1, g0 = g0, w1 = m0 / (1 - mu1), w0 = m0 / (1 - mu0)))
  
}

################################################################################

## Lasso Regression

Lasso <- function(dat) {
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV)
  
  # Step 1: Randomly split the dataset into two equal parts 
  # Double CrossFit
  
  s <- sample(1:nrow(TNDdat), nrow(TNDdat) / 2)
  TNDdat_train1 <- TNDdat[s, ]
  TNDdat_train2 <- TNDdat[-s, ]
  
  # Step 2: Estimate the probability P_TND(V = v/ C = c, Y = 0)  
  ## Training the Lasso model
  
  ### First training set: TNDdat_ctr1 
  
  TNDdat_ctr1 <- subset(TNDdat_train1, Y == 0)
  
  mod_g1_ctr <- glmnet( # The cv.glmnet function requires a set of lambda values 
    #to perform cross-validation in order to select the optimal value 
    
    y = TNDdat_ctr1$V,
    x = as.matrix(subset(TNDdat_ctr1, select = -V)),
    alpha = 1, # Lasso regularization
    lambda = cv.glmnet(y = TNDdat_ctr1$V,
                       x = as.matrix(subset(TNDdat_ctr1, select = -V)),
                       alpha = 1,
                       nfolds = 5,
                       penality.factor = 0,
                       family = "binomial")$lambda.min, # Cross-validation
    
    family = "binomial", # Categorical dependent variable
    penality.factor = 0
  )
  
  ### Second training set : TNDdat_ctr2
  
  TNDdat_ctr2 <- subset(TNDdat_train2, Y == 0)
  
  mod_g2_ctr <- glmnet(
    
    y = TNDdat_ctr2$V,
    x = as.matrix(subset(TNDdat_ctr2, select = -V)),
    alpha = 1, 
    lambda = cv.glmnet(y = TNDdat_ctr2$V,
                       x = as.matrix(subset(TNDdat_ctr2, select = -V)),
                       alpha = 1,
                       penality.factor = 0,
                       family = "binomial")$lambda.min, # Cross-validation,
    
    family = "binomial",
    penality.factor = 0
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  g1_cont <- rep(NA, nrow(TNDdat))  
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_test1 <- data.matrix(as.data.frame(cbind(select(TNDdat_train2, !c(V))) 
  ))
  
  
  TNDdata_test2 <- data.matrix(as.data.frame(cbind(select(TNDdat_train1, !c(V)) 
  )))
  
  # Predict on TNDdata_test1 (a dataset other than the one used for the model's initial training)
  
  g1_cont[-s] <- predict(mod_g1_ctr, newx = TNDdata_test1, type = "response")
  
  # Predict on TNDdata_test2 (a dataset other than the one used for the second model training)
  
  g1_cont[s] <- predict(mod_g2_ctr, newx = TNDdata_test2, type = "response")
  
  # Step 3: Estimate the probability  P_TND(Y = 1/ V = v, C = c)  
  ## Training the random forest model
  
  ### First training set
  
  Out_mu1 <- glmnet(
    
    y = TNDdat_train1$Y,
    x = data.matrix(subset(TNDdat_train1, select = -Y)),
    alpha = 1, 
    lambda = cv.glmnet(y = TNDdat_train1$Y,
                       x = data.matrix(subset(TNDdat_train1, select = -Y)),
                       alpha = 1,
                       family = "binomial")$lambda.min, # Cross-validation,
    
    family = "binomial" 
    
  )
  
  ### Second training set
  
  Out_mu2 <- glmnet(
    
    y = TNDdat_train2$Y,
    x = data.matrix(subset(TNDdat_train2, select = -Y)),
    alpha = 1, 
    lambda = cv.glmnet(y = TNDdat_train2$Y,
                       x = data.matrix(subset(TNDdat_train2, select = -Y)),
                       alpha = 1,
                       family = "binomial")$lambda.min, # Cross-validation,
    
    family = "binomial"
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  mu1 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 1, C = c) 
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_mu1_test1 <- data.matrix(as.data.frame(cbind(V = 1, select(TNDdat_train2, !c(Y, V)))))
  
  TNDdata_mu1_test2 <- data.matrix(as.data.frame(cbind(V = 1, select(TNDdat_train1, !c(Y, V)))))
  
  # Predict mu1: P(Y = 1/ V = 1) on TNDdata_mu1_test1
  
  mu1[-s] <- predict(Out_mu1, newx = TNDdata_mu1_test1, type = "response")
  
  # Predict mu1: P(Y = 1/ V = 1) on TNDdata_mu1_test2
  
  mu1[s] <- predict(Out_mu2, newx = TNDdata_mu1_test2, type = "response")
  
  ## Predicting probabilities on test sets
  # Store results
  
  mu0 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 0, C = c) 
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_mu1_test1 <- data.matrix(as.data.frame(cbind(V = 0, select(TNDdat_train2, !c(Y, V)))))
  
  TNDdata_mu1_test2 <- data.matrix(as.data.frame(cbind(V = 0, select(TNDdat_train1, !c(Y, V)))))
  
  # Predict mu0: P(Y = 1/ V = 0) on TNDdata_mu1_test1
  
  mu0[-s] <- predict(Out_mu1, newx = TNDdata_mu1_test1, type = "response")
  
  # Predict mu0: P(Y = 1/ V = 0) on TNDdata_mu1_test2
  
  mu0[s] <- predict(Out_mu2, newx = TNDdata_mu1_test2, type = "response")
  
  # Step 4: Estimate the probability  m0 (1 - Y or P(Y = 0)) 
  ## Training the Lasso model
  
  ### First training set  
  
  Out_m1 <- glmnet(
    
    y = TNDdat_train1$Y,
    x = cbind(data.matrix(subset(TNDdat_train1, select = -c(Y, V))), 0),
    alpha = 1, 
    lambda = cv.glmnet(y = TNDdat_train1$Y,
                       x = cbind(data.matrix(subset(TNDdat_train1, select = -c(Y, V))), 0),
                       alpha = 1,
                       penality.factor = 0,
                       family = "binomial")$lambda.min, # Cross-validation,
    
    family = "binomial",
    penality.factor = 0
    
  )
  
  ### Second training set
  
  Out_m2 <- glmnet(
    
    y = TNDdat_train2$Y,
    x = cbind(data.matrix(subset(TNDdat_train1, select = -c(Y, V))), 0),
    alpha = 1, 
    lambda = cv.glmnet(y = TNDdat_train2$Y,
                       x = cbind(data.matrix(subset(TNDdat_train1, select = -c(Y, V))), 0),
                       alpha = 1,
                       penality.factor = 0,
                       family = "binomial")$lambda.min, # Cross-validation,
    
    family = "binomial",
    penality.factor = 0
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  m0 <- rep(NA, nrow(TNDdat))
  
  m0[-s] <- 1 - predict(Out_m1, newx = cbind(data.matrix(select(TNDdat_train2, !c(V,Y))), 0), type = "response")
  m0[s] <- 1 - predict(Out_m2, newx = cbind(data.matrix(select(TNDdat_train1, !c(V,Y))), 0), type = "response")
  
  mu1 <- pmin(pmax(mu1, 0.01), 0.99)
  mu0 <- pmin(pmax(mu0, 0.01), 0.99)
  m0 <- pmin(pmax(m0, 0.01), 0.99)
  g1 <- pmin(pmax(g1_cont, 0.01), 0.99)
  g0 <- 1 - pmin(pmax(g1_cont, 0.01), 0.99)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1, g0 = g0, w1 = m0 / (1 - mu1), w0 = m0 / (1 - mu0)))
  
}

################################################################################

## Multivariate Adaptive Regression Splines (MARS)

Mars <- function(dat){
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV)
  
  # Step 1: Randomly split the dataset into two equal parts 
  # Double CrossFit
  
  s <- sample(1:nrow(TNDdat), nrow(TNDdat) / 2)
  
  TNDdat_train1 <- TNDdat[s, ]
  TNDdat_train2 <- TNDdat[-s, ]
  
  # Step 2: Estimate the probabilities P_TND(V = v, C = c, Y = 0)  
  ## Training the MARS model
  
  ### First training set
  
  TNDdat_train_ctr1 <- subset(TNDdat_train1, Y == 0)
  
  mod_g1_ctr <- earth(
    
    V ~ .,
    data = subset(TNDdat_train_ctr1, select = -Y),
    glm = list(family = binomial)
    
  )
  
  ### Second training set
  
  TNDdat_train_ctr2 <- subset(TNDdat_train2, Y == 0)
  
  mod_g2_ctr <- earth(
    
    V ~ .,
    data = subset(TNDdat_train_ctr2, select = -Y),
    glm = list(family = binomial)
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  g1_cont <- rep(NA, nrow(TNDdat)) 
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_test1 <- as.data.frame(cbind(select(TNDdat_train2, !c(V,Y)),
                                       V = rep(1, nrow(TNDdat_train2) ), 
                                       Y = TNDdat_train2$Y))
  
  TNDdata_test2 <- as.data.frame(cbind(select(TNDdat_train1, !c(V,Y)),
                                       V = rep(1, nrow(TNDdat_train1) ), 
                                       Y = TNDdat_train1$Y))
  
  # Predict on TNDdata_test1 (a dataset other than the one used for the model's initial training)
  
  g1_cont[-s] <- predict(mod_g1_ctr, type = "response", newdata = TNDdata_test1)
  
  # Predict on TNDdata_test2 (a dataset other than the one used for the second model training)
  
  g1_cont[s] <- predict(mod_g2_ctr, type = "response", newdata = TNDdata_test2)
  
  # Step 3: Estimate the probabilities P_TND(Y = 1/ V = v, C = c)  
  ## Training the MARS model
  
  ### First training set
  
  Out_mu1 <- earth(
    
    Y ~ .,
    data = TNDdat_train1,
    glm = list(family = binomial),
    degree = 2
    
  )
  
  ### Second training set
  
  Out_mu2 <- earth(
    
    Y ~ .,
    data = TNDdat_train2,
    glm = list(family = binomial),
    degree = 2
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  mu1 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 1, C = c)
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_mu1_test1 <- as.data.frame(cbind(V = 1, select(TNDdat_train2, !c(V,Y)) ))
  
  TNDdata_mu1_test2 <- as.data.frame(cbind(V = 1, select(TNDdat_train1, !c(V,Y)) ))
  
  # Predict mu1: P(Y = 1/ V = 1) on TNDdata_mu1_test1
  
  mu1[-s] <- predict(Out_mu1, newdata = TNDdata_mu1_test1, type = "response")
  
  # Predict mu1: P(Y = 1/ V = 1) on TNDdata_mu1_test2
  
  mu1[s] <- predict(Out_mu2, newdata = TNDdata_mu1_test2, type = "response")
  
  ## Predicting probabilities on test sets
  # Store results
  
  mu0 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 0, C = c)
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_mu1_test1 <- as.data.frame(cbind(V = 0, select(TNDdat_train2, !c(V,Y)) ))
  
  TNDdata_mu1_test2 <- as.data.frame(cbind(V = 0, select(TNDdat_train1, !c(V,Y)) ))
  
  # Predict mu0: P(Y = 1/ V = 0) on TNDdata_mu1_test1
  
  mu0[-s] <- predict(Out_mu1, newdata = TNDdata_mu1_test1, type = "response")
  
  # Predict mu0: P(Y = 1/ V = 0) on TNDdata_mu1_test2
  
  mu0[s] <- predict(Out_mu2, newdata = TNDdata_mu1_test2, type = "response")
  
  # Step 4: Estimate the probabilities  m0 (1 - Y ou P(Y = 0))   
  ## Training the MARS model
  
  ### First training set
  
  Out_m1 <- earth(
    
    Y ~ .,
    data = subset(TNDdat_train1, select = -V),
    glm = list(family = binomial)
    
  )
  
  ### Second training set
  
  Out_m2 <- earth(
    
    Y ~ .,
    data = subset(TNDdat_train2, select = -V),
    glm = list(family = binomial)
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  m0 <- rep(NA, nrow(TNDdat))
  
  m0[-s] <- 1 - predict(Out_m1, newdata = select(TNDdat_train2, !c(V,Y)), type = "response")
  m0[s] <- 1 - predict(Out_m2, newdata = select(TNDdat_train1, !c(V,Y)), type = "response")
  
  mu1 <- pmin(pmax(mu1, 0.01), 0.99)
  mu0 <- pmin(pmax(mu0, 0.01), 0.99)
  m0 <- pmin(pmax(m0, 0.01), 0.99)
  g1 <- pmin(pmax(g1_cont, 0.01), 0.99)
  g0 <- 1 - pmin(pmax(g1_cont, 0.01), 0.99)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1,g0 = g0, w1 = m0 / (1 - mu1),w0 = m0 / (1 - mu0)))
  
}

################################################################################

## Neural networks

RN <- function(dat) {
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV)
  
  # Convert the variables V and Y (response variables) into factors for classification
  
  TNDdat$V <- as.factor(TNDdat$V)
  TNDdat$Y <- as.factor(TNDdat$Y)
  
  # Step 1: Randomly split the dataset into two equal parts 
  # Double CrossFit
  
  s <- sample(1:nrow(TNDdat), nrow(TNDdat) / 2)
  TNDdat_train1 <- TNDdat[s, ]
  TNDdat_train2 <- TNDdat[-s, ]
  
  # Step 2: Estimate the probabilities P_TND(V = v/ C = c, Y = 0)  
  ## Training the RN model
  
  ### First training set: TNDdat_ctr1 
  
  TNDdat_ctr1 <- subset(TNDdat_train1, Y == 0) # Subset : controls (Y == 0)
  
  mod_g1_ctr <- nnet( # The default activation function is the sigmoid function
    # his helps the network introduce nonlinearity
    V ~ .,
    data = subset(TNDdat_ctr1, select = -Y),
    size = 5, # The number of nodes in the hidden layer
    maxit = 50, # This parameter sets the maximum number of iterations for training
    trace = FALSE
    
  )
  
  #### Second training set: TNDdat_ctr2
  
  TNDdat_ctr2 <- subset(TNDdat_train2, Y == 0) # Subset : controls (Y == 0)
  
  mod_g2_ctr <- nnet(
    
    V ~ .,
    data = subset(TNDdat_ctr2, select = -Y),
    size = 5, 
    maxit = 50,
    trace = FALSE
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  g1_cont <- rep(NA, nrow(TNDdat))  
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_test1 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_test2 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Predict on TNDdata_test1 (a dataset other than the one used for the model's initial training)
  
  g1_cont[-s] <- predict(mod_g1_ctr, newdata = TNDdata_test1, type = "raw")
  
  # Predict on TNDdata_test2 (a dataset other than the one used for the second model training)
  
  g1_cont[s] <- predict(mod_g2_ctr, newdata = TNDdata_test2, type = "raw")
  
  # Step 3: Estimate the probabilities  P_TND(Y = 1/ V = v, C = c)  
  ## Training the RN model
  
  ### First training set
  
  Out_mu1 <- nnet(
    
    Y ~ .,
    data = TNDdat_train1,
    size = 5, 
    maxit = 50,
    trace = FALSE
    
  )
  
  ### Second training set
  
  Out_mu2 <- nnet(
    
    Y ~ .,
    data = TNDdat_train2,
    size = 5, maxit = 50,
    trace = FALSE
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  mu1 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 1, C = c) 
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_mu1_test1 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_mu1_test2 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Predict mu1: P(Y = 1/ V = 1) on TNDdata_mu1_test1
  
  mu1[-s] <- predict(Out_mu1, newdata = TNDdata_mu1_test1, type = "raw")
  
  # Predict mu1: P(Y = 1/ V = 1) on TNDdata_mu1_test2
  
  mu1[s] <- predict(Out_mu2, newdata = TNDdata_mu1_test2, type = "raw")
  
  ## Predicting probabilities on test sets
  # Store results
  
  mu0 <- rep(NA, nrow(TNDdat)) #P_TND(Y = 1/ V = 0, C = c)
  
  # Make sure that the training and test sets have the same structure every time
  
  TNDdata_mu0_test1 <- TNDdat_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat_train2)), levels = levels(TNDdat$V)))
  
  TNDdata_mu0_test2 <- TNDdat_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdat_train1)), levels = levels(TNDdat$V)))
  
  # Predict mu0: P(Y = 1/ V = 0) on TNDdata_mu0_test1
  
  mu0[-s] <- predict(Out_mu1, newdata = TNDdata_mu0_test1, type = "raw")
  
  # Predict mu0: P(Y = 1/ V = 0) on TNDdata_mu0_test2
  
  mu0[s] <- predict(Out_mu2, newdata = TNDdata_mu0_test2, type = "raw")
  
  # Step 4: Estimeate the probabilities  m0 (1 - Y ou P(Y = 0))  
  ## Training the RN model
  
  ### First training set   
  
  Out_m1 <- nnet(
    
    Y ~ .,
    data = subset(TNDdat_train1, select = -V),
    size = 5,
    maxit = 50,
    trace = FALSE
    
  )
  
  ### Second training set
  
  Out_m2 <- nnet(
    
    Y ~ .,
    data = subset(TNDdat_train2, select = -V),
    size = 5, 
    maxit = 50,
    trace = FALSE
    
  )
  
  ## Predicting probabilities on test sets
  # Store results
  
  m0 <- rep(NA, nrow(TNDdat))
  
  m0[-s] <- 1 - predict(Out_m1, newdata = select(TNDdat_train2, -c(V, Y)), type = "raw")
  m0[s] <- 1 - predict(Out_m2, newdata = select(TNDdat_train1, -c(V, Y)), type = "raw")
  
  mu1 <- pmin(pmax(mu1, 0.01), 0.99)
  mu0 <- pmin(pmax(mu0, 0.01), 0.99)
  m0 <- pmin(pmax(m0, 0.01), 0.99)
  g1 <- pmin(pmax(g1_cont, 0.01), 0.99)
  g0 <- 1 - pmin(pmax(g1_cont, 0.01), 0.99)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1, g0 = g0, w1 = m0 / (1 - mu1), w0 = m0 / (1 - mu0)))
  
  
}

################################################################################
################################################################################

# Define the function for the TNDDR

TNDDR <- function(dat, methode){
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV)
  estimations <- methode(dat)
  
  # Estimation of 𝜓𝑣: 𝜓𝑣 = 𝜓𝑣(ℙTND)=𝔼TND[𝜇𝑣(𝒄)*𝜔𝑣(𝒄)] with
  #  𝜇𝑣(𝒄) = ℙTND (𝑌 = 1|𝑉 = 1,𝑪 = 𝒄),
  #  𝜔𝑣(𝒄) = 𝜋𝑣(𝒄) / 𝜋0𝑣(𝒄)
  
  A.1 <- ((1 - TNDdat$Y)*(TNDdat$V - estimations$g1))/(estimations$g1* (1 -estimations$mu1))
  psi.1 <- mean(TNDdat$Y*TNDdat$V/estimations$g1 - estimations$mu1*A.1)
  
  # Estimation of 𝜓𝑣0: 𝜓𝑣0(ℙTND)=𝔼TND[𝜇𝑣0(𝒄)*𝜔𝑣0(𝒄)] with
  #  𝜇𝑣0(𝒄) = ℙTND (𝑌 = 1|𝑉 = 0,𝑪 = 𝒄),
  
  A.0 <- ((1 - TNDdat$Y)*((1-TNDdat$V) - estimations$g0))/(estimations$g0* (1 - estimations$mu0))
  psi.0 <- mean(TNDdat$Y*(1-TNDdat$V)/estimations$g0 - estimations$mu0*A.0)
  
  # Estimation of RRM: 𝜓mRR = 𝜓𝑣∕𝜓𝑣0
  
  RRm <- pmin(pmax(psi.1/psi.0, 0.001), 0.999)
  
  # Confidence Intervals
  
  ## First method: For the TNDDR confidence interval, a logarithmic transformation can be used 
  ## to improve the accuracy of the normal approximation.
  
  ## Estimation of var(ln(𝜓eif mRR)) = var(𝔼𝕀𝔽(ln(𝜓𝑣∕𝜓𝑣0 )))
  
  log_eif_RRm <-  ((TNDdat$Y*TNDdat$V/estimations$g1 - estimations$mu1*A.1 - psi.1)/psi.1) - ((TNDdat$Y*(1-TNDdat$V)/estimations$g0 - estimations$mu0*A.0 - psi.0)/ psi.0)
  var_log_eif_RRm <-  var(log_eif_RRm)/nrow(TNDdat)
  
  ## First confidence interval for RRm
  
  IC_inf1 <- exp(log(RRm) - 1.96 * sqrt(var_log_eif_RRm) )
  IC_sup1 <- exp(log(RRm) + 1.96 * sqrt(var_log_eif_RRm) )
  
  eifpsi <- (TNDdat$Y*TNDdat$V/estimations$g1 - estimations$mu1*A.1 - psi.1)/psi.0 - RRm*(TNDdat$Y*(1-TNDdat$V)/estimations$g0 - estimations$mu0*A.0 - psi.0)/psi.0
  var <- var(eifpsi)/nrow(TNDdat)
  
  ## Second method: Second confidence interval for RRm
  
  IC_inf2 <- RRm - 1.96 * sqrt(var)
  IC_sup2 <- RRm + 1.96 * sqrt(var)
  
  ## Third method: WALD Confidence interval 
  
  varn2 <- mean((RRm* (TNDdat$Y*(1-TNDdat$V)/estimations$g0 - estimations$mu0*A.0) - (TNDdat$Y*TNDdat$V/estimations$g1 - estimations$mu1*A.1) )^2)
  denJ <- psi.0^2
  var2 <- varn2/(denJ * nrow(TNDdat))
  
  ### Third confidence interval for RRm
  
  IC_inf3 <- RRm - 1.96 * sqrt(var2)
  IC_sup3 <- RRm + 1.96 * sqrt(var2)
  
  l <- list(RRm, 1 - RRm, var_log_eif_RRm, IC_inf1, IC_sup1, 
            
            IC_inf2, IC_sup2, IC_inf3, IC_sup3)
  
  return(l)
  
}
