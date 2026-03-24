
# Install and load required packages

library(dplyr)
library(geepack)
library(geex)

# Define function for Inverse Probability Weighting (IPW)

IPW <- function(dat){
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV)
  
  # Calculate Weights 
  
  mod.denom <- glm(V ~ C,
                   
                   family = binomial(link = "logit"),  
                   data = TNDdat,
                   subset = (TNDdat$Y == 0)) # Among controls
  
  g1 <- predict(mod.denom, newdata = TNDdat, type = "response")
  
  # IPW estimator
  
  RRm <- mean(TNDdat$Y*TNDdat$V/g1)/mean(TNDdat$Y*(1-TNDdat$V)/(1-g1))
   
######## Calculate confidence intervals using the m-estimator approach #########

# Estimate variance

#### Estimation of RRm with GEEx ####

  geex_ef <- function(data){
    
    Y <- data$Y
    V <- data$V
    C <- data$C
  
    function(theta){
    
      alpha <- theta[1:2]  
    
      pscore <- (Y == 0)*plogis(alpha[1] + alpha[2]*C) # pscore only among controls
    
      # Estimation equations: weights are estimated using simple logistic regression
    
      eq_1 <- (Y == 0)*(V - pscore) # ∂l(β)/β0 = 0
    
      eq_2 <- (Y == 0)*(V - pscore)*C # ∂l(β)/β1 = 0
    
      eq_3 <- (Y*V/plogis(alpha[1] + alpha[2]*C)) - theta[3]
      eq_4 <- (Y*(1 - V)/(1 - plogis(alpha[1] + alpha[2]*C))) - theta[4]
    
      return(c(eq_1, eq_2, eq_3, eq_4))
    }
  }

  mestr <- m_estimate(estFUN = geex_ef,                                       
                    data = TNDdat,                                             
                    root_control = setup_root_control(start = c(0, 0, 0.5, 0.5))) # The same values as in the code you sent me

  beta_geex <- roots(mestr) # theta = (β0, β1, ψ10)         
  se_geex <- sqrt(diag(vcov(mestr))) # vcov(mestr) : Variance-covariance matrix
  
  # var(ln(O/E)) ~  (1/o^2)*var(O) + (1/e^2)*var(E) - (2/o*e)*COV(O,E) : O and E are dependents

  var_log_RRm <- (1/beta_geex[[3]]^2) * vcov(mestr)[3, 3] + (1/beta_geex[[4]]^2) * vcov(mestr)[4, 4] -
                
                 (2/beta_geex[[3]]*beta_geex[[4]])*vcov(mestr)[3, 4]
                

## Confidence interval for RRm

  IC_Inf <- beta_geex[[3]]/beta_geex[[4]]*exp(- 1.96*sqrt(var_log_RRm)) 
  IC_Sup <- beta_geex[[3]]/beta_geex[[4]]*exp(1.96*sqrt(var_log_RRm)) 
  
  # Store results 
  
  ## RRm (marginal relative risk)
  ## ^VE (vaccine effectiveness)
  ## Lower bound of the confidence interval
  ## Upper bound of the confidence interval

  l <- list(RRm, 1 - RRm, var_log_RRm, IC_Inf, IC_Sup)

  return(l)
   
}