# Chargement des librairies nécessaires

library(ggplot2)
library(patchwork)
library(kableExtra)

# Chargement des résultats pour les trois scénarios

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario01.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario02.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario03.RData")

# Tableau de comparaison 

comparaison01 <- data.frame(matrix(ncol = 5, 
                                   nrow = 30))

colnames(comparaison01) <- c(
  "Estimation", 
  "Performance",
  "Scénario01",
  "Scénario02",
  "Scénario03"
)

comparaison01$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison01$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison01$Scénario01 <- c(
  
  Tab1$Autres, Tab1$`Erreur de Monte Carlo`[4],
  Tab2$Autres, Tab2$`Erreur de Monte Carlo`[4],
  Tab3$Autres[1:4], Tab3$`Erreur de Monte Carlo`[4],
  Tab3$Autres[5:8], Tab3$`Erreur de Monte Carlo`[8],
  Tab3$Autres[9:12], Tab3$`Erreur de Monte Carlo`[12],
  Tab3$Autres[13:16], Tab3$`Erreur de Monte Carlo`[16]
  
)

comparaison01$Scénario02 <- c(
  
  Tab1_2$Autres, Tab1_2$`Erreur de Monte Carlo`[4],
  Tab2_2$Autres, Tab2_2$`Erreur de Monte Carlo`[4],
  Tab3_2$Autres[1:4], Tab3_2$`Erreur de Monte Carlo`[4],
  Tab3_2$Autres[5:8], Tab3_2$`Erreur de Monte Carlo`[8],
  Tab3_2$Autres[9:12], Tab3_2$`Erreur de Monte Carlo`[12],
  Tab3_2$Autres[13:16], Tab3_2$`Erreur de Monte Carlo`[16]
  
)

comparaison01$Scénario03 <- c(
  
  Tab1_3$Autres, Tab1_3$`Erreur de Monte Carlo`[4],
  Tab2_3$Autres, Tab2_3$`Erreur de Monte Carlo`[4],
  Tab3_3$Autres[1:4], Tab3_3$`Erreur de Monte Carlo`[4],
  Tab3_3$Autres[5:8], Tab3_3$`Erreur de Monte Carlo`[8],
  Tab3_3$Autres[9:12], Tab3_3$`Erreur de Monte Carlo`[12],
  Tab3_3$Autres[13:16], Tab3_3$`Erreur de Monte Carlo`[16]
  
)

# nsim = 1000, couverture vaccinale  ~ 33% : Scénarios 01, 02 et 03

comparaison01$Scénario01 <- sapply(comparaison01$Scénario01, FUN = function(x)x[[2]])
comparaison01$Scénario02 <- sapply(comparaison01$Scénario02, FUN = function(x)x[[2]][[1]])
comparaison01$Scénario03 <- sapply(comparaison01$Scénario03, FUN = function(x)x[[2]][[1]])

View(comparaison01)

# nsim = 1000, coouverture vaccinale  ~ 33%, 50% et 70% : Scénarios 02 et 03

comparaison02 <- data.frame(matrix(ncol = 8, 
                                   nrow = 30))

colnames(comparaison02) <- c(
  "Estimation", 
  "Performance",
  "Scénario02_33%",
  "Scénario02_50%",
  "Scénario02_70%",
  "Scénario03_33%",
  "Scénario03_50%",
  "Scénario03_70%"
)

comparaison02$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison02$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison02$`Scénario02_33%` <- c(
  
  sapply(Tab1_2$Autres, function(x) x[[2]][[1]]), Tab1_2$`Erreur de Monte Carlo`[[4]][[2]][[1]],
  sapply(Tab1_2$Autres, function(x) x[[2]][[1]]), Tab2_2$`Erreur de Monte Carlo`[[4]][[2]][[1]],
  sapply(Tab3_2$Autres[1:4],function(x) x[[2]][[1]]) , Tab3_2$`Erreur de Monte Carlo`[[4]][[2]][[1]],
  sapply(Tab3_2$Autres[5:8],function(x) x[[2]][[1]]), Tab3_2$`Erreur de Monte Carlo`[[8]][[2]][[1]],
  sapply(Tab3_2$Autres[9:12],function(x) x[[2]][[1]]), Tab3_2$`Erreur de Monte Carlo`[[12]][[2]][[1]],
  sapply(Tab3_2$Autres[13:16],function(x) x[[2]][[1]]), Tab3_2$`Erreur de Monte Carlo`[[16]][[2]][[1]]
  
)

comparaison02$`Scénario02_50%` <- c(
  
  sapply(Tab1_2$Autres, function(x) x[[2]][[2]]), Tab1_2$`Erreur de Monte Carlo`[[4]][[2]][[2]],
  sapply(Tab1_2$Autres, function(x) x[[2]][[2]]), Tab2_2$`Erreur de Monte Carlo`[[4]][[2]][[2]],
  sapply(Tab3_2$Autres[1:4],function(x) x[[2]][[2]]) , Tab3_2$`Erreur de Monte Carlo`[[4]][[2]][[2]],
  sapply(Tab3_2$Autres[5:8],function(x) x[[2]][[2]]), Tab3_2$`Erreur de Monte Carlo`[[8]][[2]][[2]],
  sapply(Tab3_2$Autres[9:12],function(x) x[[2]][[2]]), Tab3_2$`Erreur de Monte Carlo`[[12]][[2]][[2]],
  sapply(Tab3_2$Autres[13:16],function(x) x[[2]][[2]]), Tab3_2$`Erreur de Monte Carlo`[[16]][[2]][[2]]
  
)

comparaison02$`Scénario02_70%` <- c(
  
  sapply(Tab1_2$Autres, function(x) x[[2]][[3]]), Tab1_2$`Erreur de Monte Carlo`[[4]][[2]][[3]],
  sapply(Tab1_2$Autres, function(x) x[[2]][[3]]), Tab2_2$`Erreur de Monte Carlo`[[4]][[2]][[3]],
  sapply(Tab3_2$Autres[1:4],function(x) x[[2]][[3]]) , Tab3_2$`Erreur de Monte Carlo`[[4]][[2]][[3]],
  sapply(Tab3_2$Autres[5:8],function(x) x[[2]][[3]]), Tab3_2$`Erreur de Monte Carlo`[[8]][[2]][[3]],
  sapply(Tab3_2$Autres[9:12],function(x) x[[2]][[3]]), Tab3_2$`Erreur de Monte Carlo`[[12]][[2]][[3]],
  sapply(Tab3_2$Autres[13:16],function(x) x[[2]][[3]]), Tab3_2$`Erreur de Monte Carlo`[[16]][[2]][[3]]
  
)

comparaison02$`Scénario03_33%` <- c(
  
  sapply(Tab1_3$Autres, function(x) x[[2]][[1]]), Tab1_3$`Erreur de Monte Carlo`[[4]][[2]][[1]],
  sapply(Tab1_3$Autres, function(x) x[[2]][[1]]), Tab2_3$`Erreur de Monte Carlo`[[4]][[2]][[1]],
  sapply(Tab3_3$Autres[1:4],function(x) x[[2]][[1]]) , Tab3_3$`Erreur de Monte Carlo`[[4]][[2]][[1]],
  sapply(Tab3_3$Autres[5:8],function(x) x[[2]][[1]]), Tab3_3$`Erreur de Monte Carlo`[[8]][[2]][[1]],
  sapply(Tab3_3$Autres[9:12],function(x) x[[2]][[1]]), Tab3_3$`Erreur de Monte Carlo`[[12]][[2]][[1]],
  sapply(Tab3_3$Autres[13:16],function(x) x[[2]][[1]]), Tab3_3$`Erreur de Monte Carlo`[[16]][[2]][[1]]
  
)

comparaison02$`Scénario03_50%` <- c(
  
  sapply(Tab1_3$Autres, function(x) x[[2]][[2]]), Tab1_3$`Erreur de Monte Carlo`[[4]][[2]][[2]],
  sapply(Tab1_3$Autres, function(x) x[[2]][[2]]), Tab2_3$`Erreur de Monte Carlo`[[4]][[2]][[2]],
  sapply(Tab3_3$Autres[1:4],function(x) x[[2]][[2]]) , Tab3_3$`Erreur de Monte Carlo`[[4]][[2]][[2]],
  sapply(Tab3_3$Autres[5:8],function(x) x[[2]][[2]]), Tab3_3$`Erreur de Monte Carlo`[[8]][[2]][[2]],
  sapply(Tab3_3$Autres[9:12],function(x) x[[2]][[2]]), Tab3_3$`Erreur de Monte Carlo`[[12]][[2]][[2]],
  sapply(Tab3_3$Autres[13:16],function(x) x[[2]][[2]]), Tab3_3$`Erreur de Monte Carlo`[[16]][[2]][[2]]
  
)

comparaison02$`Scénario03_70%` <- c(
  
  sapply(Tab1_3$Autres, function(x) x[[2]][[3]]), Tab1_3$`Erreur de Monte Carlo`[[4]][[2]][[3]],
  sapply(Tab1_3$Autres, function(x) x[[2]][[3]]), Tab2_3$`Erreur de Monte Carlo`[[4]][[2]][[3]],
  sapply(Tab3_3$Autres[1:4],function(x) x[[2]][[3]]) , Tab3_3$`Erreur de Monte Carlo`[[4]][[2]][[3]],
  sapply(Tab3_3$Autres[5:8],function(x) x[[2]][[3]]), Tab3_3$`Erreur de Monte Carlo`[[8]][[2]][[3]],
  sapply(Tab3_3$Autres[9:12],function(x) x[[2]][[3]]), Tab3_3$`Erreur de Monte Carlo`[[12]][[2]][[3]],
  sapply(Tab3_3$Autres[13:16],function(x) x[[2]][[3]]), Tab3_3$`Erreur de Monte Carlo`[[16]][[2]][[3]]
  
)

View(comparaison02)
