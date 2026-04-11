
# Install and load required packages

library("stats")
library("rje")


datagen <- function(seed = sample(1:1000000, size = 1), popsize = 1*10**6, ssize = 5000, 
                    
                    co_inf_para1 = 0, co_inf_para2 = 0,
                    
                    CV = 0.33, # Vaccination coverage
                    I1_prev = 0.15, # Prevalence of I1
                    I2_prev = 0.5, # Prevalence of I2
                    W1_prev = 0.05,# Prevalence of W1
                    W2_prev = 0.02, # Prevalence of W2
                                   # Default values for the baseline scenario
                    
                    return_full = FALSE) { # To choose between the population and the TND sample
  
  
  # Generation of the continuous confounding factor C 
  
  C <- runif(n = popsize, 0.1, 3) # We start with the root nodes of our DAG
  
  # Generate Vaccination Status V : V ~ Bernoulli(logit(b0 + b1*C))
  
  ## expressing the desired marginal probability as a function of b0 for V
  
  a <- 1/0.87
  f_V <- function (b0) {
    
    a*(log(1 + exp(b0 + 3*0.3)) - log(1 + exp(b0 + 0.1*0.3))) - CV
      
       # CV: the desired marginal probability, Vaccination coverage (hyperparameter)
    
  }
  
  b0_V <- uniroot(f_V, c(-10, 10))$root
  
  V <- rbinom(n = popsize, size = 1, prob = plogis(b0_V + 0.3*C)) 
  # See the file Anciens fichiers\Génération des variables
  # for the value of b0 
  
  # Generate independent infections I1 : I1 ~ Bernoulli(logit(b0 + b1*C + b2*C2)) and 
  #                                 I2 : I2 ~ Bernoulli(logit(b0 + b1*C + b2*V + b3*C2))
  
  C2 <- 2*rbinom(n = popsize, 1, 0.5) - 1
  
  ## expressing the desired marginal probability as a function of b0 for I1 with b1 and b2 known
  
  a <- 0.5/2.9
  
  f_I1 <- function (b0) { # P(I1 = 1) - 0.15
    
    (a/0.35)*(log(1 + exp(b0 + 3*0.35 - co_inf_para1)) + log(1 + exp(b0 + 3*0.35 + co_inf_para1)) - 
              
    log(1 + exp(b0 + 0.1*0.35 - co_inf_para1)) - log(1 + exp(b0 + 0.1*0.35 + co_inf_para1))) - 
      
    I1_prev # I1_prev: the desired marginal probability, prevalence of I1 (hyperparameter)
    
  }
  
  b0 <- uniroot(f_I1, c(-10, 10))$root
  
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.35*C + co_inf_para1*C2)) 
  
  ## expressing the desired marginal probability as a function of b0 for I2 with b1, b2 and b3 known
  
  P_I2 <- function(c,b0) { # Joint probability for I2
    
    a*(expit(b0 + 0.15*c - 0.1 - co_inf_para2)*expit(b0_V + 0.3*c) + 
         
         expit(b0 + 0.15*c - co_inf_para2)*(1 - expit(b0_V + 0.3*c)) + 
         
         expit(b0 + 0.15*c - 0.1 + co_inf_para2)*expit(b0_V + 0.3*c) + 
         
         expit(b0 + 0.15*c + co_inf_para2)*(1 - expit(b0_V + 0.3*c))
    )
    
  }
  
  f_I2 <- function(b0){ # P(I2 = 2) - 0.50
    
    integrate(P_I2, 0.1, 3, b0)$value - I2_prev # I2_prev: the desired marginal probability, prevalence of I2 (hyperparameter)
    
  }
  
  b0 <- uniroot(f_I2, c(-10, 10))$root
  
  I2 <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.15*C - 0.1*V + co_inf_para2*C2)) # Pour une prévalence ~ 
  
  # Generate symptoms W1: W1~Bernoulli(logit(b0 + b1*C[I1 = 1]))  and 
  #                   W2: W2~Bernoulli(logit(b0 + b1*C[I2 == 1] + b2*V[I2 == 1]))
  
  ## expressing the desired marginal probability as a function of b0 for W1
  
  a <- 1/2.9
  b1 <- 0.5
  
  f_W1 <- function (b0) {
    
    (a/b1)*(log(1 + exp(b0 + 3*b1 ))  - log(1 + exp(b0 + 0.1*b1))) - 
      
      W1_prev # W1_prev: the desired marginal probability, prevalence of W1 (hyperparameter)
    
  }
  
  b0 <- uniroot(f_W1, c(-10, 10))$root 
  
  W1 <- rep(0, popsize) # W1 = 0 if I1 = 0, so:
  
  W1[I1 == 1] <- rbinom(
    
    n = sum(I1 == 1),
    size = 1,
    prob = plogis(b0 + 0.5 * C[I1 == 1])
    
  ) # W1 = 1 ~ 4%-5%
  
  ## expressing the desired marginal probability as a function of b0 for W2
  
  a <- 0.5/2.9

  b1 <- 2
  b2 <- -0.91
  
  P_W2 <- function(c,b0) { # Joint probability for W2
    
    a*( expit(b0 + b1*c + b2)*expit(b0_V + 0.3*c) + 
         
         expit(b0 + b1*c)*(1 - expit(b0_V + 0.3*c))
    )
    
  }
  
  help(integrate)
  # ...	
  # Additional arguments to be passed to f: this argument is used to pass 
  # the value b0 to the P_W2 function in this case.
  
  f_W2 <- function(b0){
    
    integrate(P_W2, 0.1, 3, b0)$value - W2_prev # W2_prev: the desired marginal probability, prevalence of W2 (hyperparameter)
    
  }
  
  b0 <- uniroot(f_W2, c(-10, 10))$root 
   
  W2 <- rep(0, popsize) # W2_0 = 0 if I2_0 = 0, so:
  
  W2[I2 == 1] <- rbinom(
    
    n = sum(I2 == 1),
    size = 1,
    prob = plogis(b0 + 2*C[I2 == 1] - 0.91*V[I2 == 1])
    
  ) # W2 = 1 ~ 3%
  
  co_W <- sum(W1 == 1 & W2 == 1, na.rm = TRUE)
  per_co_W <- co_W / popsize * 100
  
  # Generate W
  
  W <- pmax(W1, W2)
  
  # Generate hospitalization
  
  H = rep(0, popsize) # H = 0 if W = 0 (only possible if symptoms present), so:
  
  H[W == 1] <- rbinom(prob = plogis(-1.5 + 0.5*C[W == 1]),
                      size = 1, n = sum(W == 1))
  
  # The TND includes only those who have been tested : H = 1
  
  R <- sample(which(W == 1), ssize, replace = TRUE) # A random sample of patients with symptoms (W=1 <=> H=1) with replacement

  if (return_full == FALSE) {
    
    dat <- as.data.frame(cbind(Infec_RSV = I2, Infec = I1, H = H, W1 = W1, W2 = W2,
                               
                               W = W, V = V, C = C)[R, ]) # Respiratory syncytial virus RSV
    
    pourcentage_unique <- (nrow(unique(dat)) / nrow(dat)) * 100 # Pourcentages des non doublons
    
    # Calculate the percentage of (symptomatic) co-infections in the sample
    
    co_inf <- sum(dat$Infec == 1 & dat$Infec_RSV == 1)
    
    per_co_inf <- co_inf / ssize * 100
    
    print(paste("The percentage of co-infections in the sample is :", per_co_inf))
    print(paste("The percentage of unique lines is :", pourcentage_unique))
    
  } else { # Data for the total population (return_full == TRUE)
    
    dat <- as.data.frame(cbind(Infec_RSV = I2, Infec = I1, H = H, W1 = W1, W2 = W2,
                               
                               W = W, V = V, C = C))
    
    # Calculate the percentage of symptomatic co-infections in the population
    
    co_inf_1 <- sum(dat$Infec == 1 & dat$Infec_RSV == 1 & dat$W == 1)
    
    per_co_inf_1 <- co_inf_1 / popsize * 100
    
    print((paste("The percentage of symptomatic co-infections in the population is :", per_co_inf_1)))
    
    # Calculate the percentage of asymptomatic co-infections in the population
    
    co_inf_0 <- sum(dat$Infec == 1 & dat$Infec_RSV == 1 & dat$W == 0)
    
    per_co_inf_0 <- co_inf_0 / popsize * 100
    
    print((paste("The percentage of asymptomatic co-infections in the population is :", per_co_inf_0)))
    
  }
  
  return(dat)
  
}  

dat_full <- datagen(return_full = TRUE)

sum(dat_full$V == 1)/nrow(dat_full) * 100 # Vaccination coverage

sum(dat_full$Infec == 1)/nrow(dat_full) * 100
sum(dat_full$W1 == 1)/sum(dat_full$Infec == 1) * 100 # Symptomatic prevalence for I1

sum(dat_full$Infec_RSV == 1)/nrow(dat_full) * 100
sum(dat_full$W2 == 1)/sum(dat_full$Infec_RSV == 1) * 100 # Symptomatic prevalence for I2

################################################################################ 

# Generating counterfactual scenarios

datagen.cont <- function(seed = sample(1:1000000, size = 1), popsize = 1*10**6,
                         
                         co_inf_para1 = 0, co_inf_para2 = 0,
                         
                         CV = 0.33, # Vaccination coverage
                         I1_prev = 0.15, # Prevalence of I1
                         I2_prev = 0.5, # Prevalence of I2
                         W1_prev = 0.05,# Prevalence of W1
                         W2_prev = 0.02, # Prevalence of W2
                         # Default values for the baseline scenario
                         
                         return_full = FALSE) { # To choose between the population and the TND sample
  
  
  # Generation of the continuous confounding factor C 
  
  C <- runif(n = popsize, 0.1, 3); 
  
  # Generate Vaccination Status V : V ~ Bernoulli(logit(b0 + b1*C))
  
  ## expressing the desired marginal probability as a function of b0 for V
  
  a <- 1/0.87
  f_V <- function (b0) {
    
    a*(log(1 + exp(b0 + 3*0.3)) - log(1 + exp(b0 + 0.1*0.3))) - 
      
      CV # CV: the desired marginal probability, Vaccination coverage (hyperparameter)
    
  }
  
  b0_V <- uniroot(f_V, c(-10, 10))$root
  
  V <- rbinom(n = popsize, size = 1, prob = plogis(b0_V + 0.3*C)) 
  # See the file Anciens fichiers\Génération des variables
  # for the value of b0 
  
  # Generate independent infections I1 : I1 ~ Bernoulli(logit(b0 + b1*C + b2*C2)) and 
  #                                 I2 : I2 ~ Bernoulli(logit(b0 + b1*C + b2*V + b3*C2))
  
  
  C2 <- 2*rbinom(n = popsize, 1, 0.5) - 1
  
  ## expressing the desired marginal probability as a function of b0 for I1 with b1 and b2 known
  
  a <- 0.5/2.9
  
  f_I1 <- function (b0) { # P(I1 = 1) - 0.15
    
    (a/0.35)*(log(1 + exp(b0 + 3*0.35 - co_inf_para1)) + log(1 + exp(b0 + 3*0.35 + co_inf_para1)) - 
                
              log(1 + exp(b0 + 0.1*0.35 - co_inf_para1)) - log(1 + exp(b0 + 0.1*0.35 + co_inf_para1))) - 
      
      I1_prev # I1_prev: the desired marginal probability, prevalence of I1 (hyperparameter)
    
  }
  
  b0 <- uniroot(f_I1, c(-10, 10))$root
  
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.35*C + co_inf_para1*C2))
  
  ## expressing the desired marginal probability as a function of b0 for I2 with b1, b2 and b3 known
  
  P_I2 <- function(c,b0) { # Joint probability for I2
    
    a*(expit(b0 + 0.15*c - 0.1 - co_inf_para2)*expit(b0_V + 0.3*c) + 
         
         expit(b0 + 0.15*c - co_inf_para2)*(1 - expit(b0_V + 0.3*c)) + 
         
         expit(b0 + 0.15*c - 0.1 + co_inf_para2)*expit(b0_V + 0.3*c) + 
         
         expit(b0 + 0.15*c + co_inf_para2)*(1 - expit(b0_V + 0.3*c))
    )
    
  }
  
  f_I2 <- function(b0){ # P(I2 = 1) - 0.50
    
    integrate(P_I2, 0.1, 3, b0)$value - I2_prev # I2_prev: the desired marginal probability, prevalence of I2 (hyperparameter)
    
  }
  
  b0 <- uniroot(f_I2, c(-10, 10))$root
  
  I2_0 <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.15*C - 0.1*0 + co_inf_para2*C2))
  I2_1 <- rbinom(n = popsize, size = 1, prob = plogis(b0 + 0.15*C - 0.1*1 + co_inf_para2*C2))
  
  # Generate symptoms W1: W1~Bernoulli(logit(b0 + b1*C[I1 = 1]))  and 
  #                   W2: W2~Bernoulli(logit(b0 + b1*C[I2 == 1] + b2*V[I2 == 1]))
  
  ## expressing the desired marginal probability as a function of b0 for W1
  
  a <- 1/2.9
  b1 <- 0.5
  
  f_W1 <- function (b0) {
    
    (a/b1)*(log(1 + exp(b0 + 3*b1 ))  - log(1 + exp(b0 + 0.1*b1))) - 
      
      W1_prev #  W1_prev: the desired marginal probability, prevalence of W1 (hyperparameter)
    
  }
  
  b0 <- uniroot(f_W1, c(-10, 10))$root
  W1 <- rep(0, popsize) # W1 = 0 if I1 = 0, so:
  
  W1[I1 == 1] <- rbinom(
    
    n = sum(I1 == 1),
    size = 1,
    prob = plogis(b0 + 0.5 * C[I1 == 1])
    
  )
  
  ## expressing the desired marginal probability as a function of b0 for W1
  
  a <- 0.5/2.9
  
  b1 <- 2
  b2 <- -0.91
  
  P_W2 <- function(c,b0) { # Joint probability for W2
    
    a*(expit(b0 + b1*c + b2)*expit(b0_V + 0.3*c) + 
         
         expit(b0 + b1*c)*(1 - expit(b0_V + 0.3*c)) 
    )
    
  }
  
  help(integrate)
  # ...	
  # Additional arguments to be passed to f: this argument is used to pass 
  # the value b0 to the P_W2 function in this case.
  
  f_W2 <- function(b0){
    
    integrate(P_W2, 0.1, 3, b0)$value - W2_prev # W2_prev: the desired marginal probability, prevalence of W2 (hyperparameter)
    
  }
  
  b0 <- uniroot(f_W2, c(-10, 10))$root
  
  W2_0 <- rep(0, popsize) # W2_0 = 0 if I2_0 = 0, so:
  W2_1 <- rep(0, popsize) # W2_1 = 0 if I2_1 = 0, soc:
  
  W2_0[I2_0 == 1] <- rbinom(
    
    n = sum(I2_0 == 1),
    size = 1,
    prob = plogis(b0 + 2*C[I2_0 == 1] - 0.91*0) # V[I2_0 == 1]
    
  )
  
  W2_1[I2_1 == 1] <- rbinom(
    
    n = sum(I2_1 == 1), 
    size = 1,
    prob = plogis(b0 + 2*C[I2_1 == 1] - 0.91*1) # V[I2_1 == 1]
    
  )
  
  # Generate W
  
  W_0 <- pmax(W1, W2_0)
  W_1 <- pmax(W1, W2_1)
  
  
  # Generate hospitalization 
  
  H_0 = rep(0, popsize) # H_0 = 0 if W_0 = 0, so:
  H_1 = rep(0, popsize) # H_1 = 0 if W_1 = 0, so:
  
  H_0[W_0 == 1] <- rbinom(prob = plogis(-1.5 + 0.5*C[W_0 == 1]),
                          size = 1, n = sum(W_0 == 1))
  
  H_1[W_1 == 1] <- rbinom(prob = plogis(-1.5 + 0.5*C[W_1 == 1]),
                          size = 1, n = sum(W_1 == 1))
  
  dat <- data.frame(C, I1, I2_0, I2_1, W1, W2_0, W2_1, W_0, W_1, H_0, H_1)
  
  # Calculate the percentage of co-infections
  
  co_inf_0 <- sum(dat$I1 == 1 & dat$I2_0 == 1)
  
  per_co_inf_0 <- co_inf_0 / popsize * 100
  
  co_inf_1 <- sum(dat$I1 == 1 & dat$I2_1 == 1)
  
  per_co_inf_1 <- co_inf_1 / popsize * 100
  
  print(paste(per_co_inf_0, per_co_inf_1))
  
  return(dat)
}
