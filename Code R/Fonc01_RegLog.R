
# Install and load required packages

library(dplyr)

# Define function for logistic regression

RegLog <- function(dat){
  
  TNDdat <- data.frame(C = dat$C, V = dat$V, Y = dat$Infec_RSV)
 
  # Analysis using logistic regression model
  
  fit.TND <- glm(Y ~ V + C,
                 
                 family = binomial(link = "logit"),  
                 # For TND data, we use the logit model
                 data = TNDdat) 
  
  # Store results
  
  resultats.TND <- summary(fit.TND)
  
  ## Logistic regression coefficient
  ## Standard error of the coefficient
  ## RRc (conditional risk ratio)
  ## ^VE (vaccine effectiveness)
  
  l <- list(resultats.TND$coefficients[2,1], resultats.TND$coefficients[2,2],
            exp(resultats.TND$coefficients[2,1]), 1 - exp(resultats.TND$coefficients[2,1]))  
  

  return(l)
  
}

