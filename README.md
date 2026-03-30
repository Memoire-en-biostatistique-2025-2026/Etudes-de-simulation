# Etudes-de-simulation

## Introduction

 This project aim to explore the impact of the proportion of co-infections on the estimation of vaccine effectiveness in the context of a testnegative design study (TND). And this by trying to isolate the impact of co-infections while maintaining fixed baseline prevalences.
 
 The simulation study aims to evaluate the performance of different methods for estimating vaccine effectiveness : logistic regression, IPW (Inverse Probability Weighting), and TNDDR (doubly robust estimator), in the presence of co-infection, which may violate the control exchangeability assumption.
 
 Next we tried to explore individual impact of different hyperparameters, including vaccination coverage, prevalence of target infection, prevalence of non-target infection and sample size.
 
## Notes

•	The intercept values for generating the study variables were calculated by expressing the desired marginal probabilities as a function of b0 for a fixed value of the other parameters. This allows the user to choose the values of the different hyperparameters: vaccination coverage, prevalence of infection, etc. while retaining fixed prevalences.  

 We want to generate P(I1 = 1| C2) = b0 + b1*C + b2C2 so that, the marginal probability P(I1 = 1) is fixed at 0.15. Because C and C2 are independent by construction, the constraint P(I1 = 1) = 0.15 can be expressed as:
 
 0.15 = ∫_0.1^3〖P(I1=1│C= c,C2 = -1) * P(C2 = -1)f(c)dc〗+ ∫_0.1^3〖P(I1= 1|C = c,C2 = 1) * P(C2 = 1)f(c) dc〗, with b1 and b2 known and b0 to be found.

•	Most values of the various hyperparameters for the base scenario were chosen based on the literature to create realistic scenarios.

•	The variable C2, which affects both I1 and I2 with two different parameters, allows co-infection to be controlled without affecting the prevalence of I1 and I2; the aim is to study the isolated effect of co-infection while maintaining the other hyperparameters fixed.

## Content

1. Source codes and data generation environment (Folder `Code R`):

   a) `Générations des données.R`: Data generation process (variables and their distributions) as well as counterfactual scenarios.
   
   b) Estimation functions:
      1. `Fonc01_RegLog.R`: Code for logistic regression.
      2. `Fonc02_IPW.R`: Code for Inverse Probability Weighting (IPW) method.
      3. `Fonc03_TNDDR.R`: Code for an efficient and doubly robust estimation of vaccine effectiveness under the TND.
      
   c) Scenarios (for n(sample size) = 1000 and n = 5000):
   
   * Scenario 01: Baseline scenario; co-infection_TND ~ 0%; co_infe_para1 = 8, co_ine_para2 = -8.

   * Scenario 02: Moderate co-infection rate_TND ~ 10%; co_infe_para1 = 2, co_ine_para2 = -2.

   * Scenario 03: High co-infection rate_TND ~ 20% co_infe_para1 = 1 and co_ine_para2 = -1.
   
   * Scenario 04: High co-infection rate_TND ~ 40% co_infe_para1 = 1 and co_infe_para2 = 2.
   
  Then we changed one hyperparameter at a time, keeping the others at their initial values (those of the baseline scenario).
  
   High co-infection rate_TND ~ 40% for:
   
   * Scenario 05: Prevalence of I2 (infection of interest): 15%, 50%, and 70%.

   * Scenario 06: Vaccination coverage: 33%, 50%, and 70%

   * Scenario 07: Prevalence of I1 (infection of interest): 10%, 30%, and 50%.
   
   * Scenario 08: Prevalence of W1 (infection of interest): 10%, 15%, and 30%.
   
   * Scenario 09: Prevalence of W2 (infection of interest): (Not coded yet)

2. `Résultats bruts`: The raw results of some of the scenarios presented above.
