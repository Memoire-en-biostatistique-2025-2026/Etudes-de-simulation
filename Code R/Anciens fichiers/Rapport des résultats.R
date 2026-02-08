# Chargement des librairies nécessaires

library(ggplot2)
library(patchwork)
library(kableExtra)

# Chargement des résultats pour les trois scénarios

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario01.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario02.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario03.RData")

# Tableau de comparaison 

comparaison <- data.frame(matrix(ncol = 5, 
                                 nrow = 30))

colnames(comparaison) <- c(
                         "Estimation", 
                         "Performance",
                         "Scénario01",
                         "Scénario02",
                         "Scénario03"
)

comparaison$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                            "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                            "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

a <- as.data.frame(c(Tab1$Autres, Tab1$`Erreur de Monte Carlo`))
b <- as.data.frame(c(Tab2$Autres, Tab2$`Erreur de Monte Carlo`))
c <- as.data.frame(c(Tab3$Autres, Tab3$`Erreur de Monte Carlo`))

comparaison$Scénario01 <- c(a$bias, a$var, a$mse, a$value, a$coverage,
                            b$bias, b$var, b$mse, b$value, b$coverage,
                            c$bias, c$var, c$mse, c$value, c$coverage,
                            c$bias.1, c$var.1, c$mse.1, c$value.1, c$coverage.1,
                            c$bias.2, c$var.2, c$mse.2, c$value.2, c$coverage.2,
                            c$bias.3, c$var.3, c$mse.3, c$value.3, c$coverage.3)
  
a <- as.data.frame(c(Tab1_2$Autres, Tab1_2$`Erreur de Monte Carlo`))
b <- as.data.frame(c(Tab2_2$Autres, Tab2_2$`Erreur de Monte Carlo`))
c <- as.data.frame(c(Tab3_2$Autres, Tab3_2$`Erreur de Monte Carlo`))

comparaison$Scénario02 <- c(a$value.bias[1], a$value.var[1], a$value.mse[1], a$value[1], a$value.coverage[1],
                            b$value.bias[1], b$value.var[1], b$value.mse[1], b$value[1], b$value.coverage[1],
                            c$value.bias[1], c$value.var[1], c$value.mse[1], c$value[1], c$value.coverage[1],
                            c$value.bias.3[1], c$value.var.3[1], c$value.mse.3[1], c$value.1[1], c$value.coverage.3[1],
                            c$value.bias.6[1], c$value.var.6[1], c$value.mse.6[1], c$value.2[1], c$value.coverage.6[1],
                            c$value.bias.9[1], c$value.var.9[1], c$value.mse.9[1], c$value.3[1], c$value.coverage.9[1])

a <- as.data.frame(c(Tab1_3$Autres, Tab1_3$`Erreur de Monte Carlo`))
b <- as.data.frame(c(Tab2_3$Autres, Tab2_3$`Erreur de Monte Carlo`))
c <- as.data.frame(c(Tab3_3$Autres, Tab3_3$`Erreur de Monte Carlo`))

comparaison$Scénario03 <- c(a$value.bias[1], a$value.var[1], a$value.mse[1], a$value[1], a$value.coverage[1],
                            b$value.bias[1], b$value.var[1], b$value.mse[1], b$value[1], b$value.coverage[1],
                            c$value.bias[1], c$value.var[1], c$value.mse[1], c$value[1], c$value.coverage[1],
                            c$value.bias.3[1], c$value.var.3[1], c$value.mse.3[1], c$value.1[1], c$value.coverage.3[1],
                            c$value.bias.6[1], c$value.var.6[1], c$value.mse.6[1], c$value.2[1], c$value.coverage.6[1],
                            c$value.bias.9[1], c$value.var.9[1], c$value.mse.9[1], c$value.3[1], c$value.coverage.9[1])
################################################################################

View(comparaison)
