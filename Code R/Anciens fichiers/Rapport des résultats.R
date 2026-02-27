# Chargement des librairies nécessaires

library(ggplot2)
library(patchwork)
library(kableExtra)
library("gridExtra")
library("cowplot")
library(ggpubr)
library(Rmisc)

# Chargement des résultats pour les trois scénarios

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario01.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario02.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario03.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario04.RData")

# Tableau de comparaison 

comparaison01 <- data.frame(matrix(ncol = 6, 
                                   nrow = 30))

colnames(comparaison01) <- c(
  "Estimation", 
  "Performance",
  "Scénario01_0%",
  "Scénario02_10%",
  "Scénario03_20%",
  "Scénario04_40%"
)

comparaison01$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison01$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison01$`Scénario01_0%` <- c(
  
  Tab_1$Autres, Tab_1$`Erreur de Monte Carlo`[4],
  Tab2_1$Autres, Tab2_1$`Erreur de Monte Carlo`[4],
  Tab3_1$Autres[1:4], Tab3_1$`Erreur de Monte Carlo`[4],
  Tab3_1$Autres[5:8], Tab3_1$`Erreur de Monte Carlo`[8],
  Tab3_1$Autres[9:12], Tab3_1$`Erreur de Monte Carlo`[12],
  Tab3_1$Autres[13:16], Tab3_1$`Erreur de Monte Carlo`[16]
  
)

comparaison01$`Scénario02_10%` <- c(
  
  Tab_2$Autres, Tab_2$`Erreur de Monte Carlo`[4],
  Tab2_2$Autres, Tab2_2$`Erreur de Monte Carlo`[4],
  Tab3_2$Autres[1:4], Tab3_2$`Erreur de Monte Carlo`[4],
  Tab3_2$Autres[5:8], Tab3_2$`Erreur de Monte Carlo`[8],
  Tab3_2$Autres[9:12], Tab3_2$`Erreur de Monte Carlo`[12],
  Tab3_2$Autres[13:16], Tab3_2$`Erreur de Monte Carlo`[16]
  
)

comparaison01$`Scénario03_20%` <- c(
  
  Tab_3$Autres, Tab_3$`Erreur de Monte Carlo`[4],
  Tab2_3$Autres, Tab2_3$`Erreur de Monte Carlo`[4],
  Tab3_3$Autres[1:4], Tab3_3$`Erreur de Monte Carlo`[4],
  Tab3_3$Autres[5:8], Tab3_3$`Erreur de Monte Carlo`[8],
  Tab3_3$Autres[9:12], Tab3_3$`Erreur de Monte Carlo`[12],
  Tab3_3$Autres[13:16], Tab3_3$`Erreur de Monte Carlo`[16]
  
)

comparaison01$`Scénario04_40%` <- c(
  
  Tab_4$Autres, Tab_4$`Erreur de Monte Carlo`[4],
  Tab2_4$Autres, Tab2_4$`Erreur de Monte Carlo`[4],
  Tab3_4$Autres[1:4], Tab3_4$`Erreur de Monte Carlo`[4],
  Tab3_4$Autres[5:8], Tab3_4$`Erreur de Monte Carlo`[8],
  Tab3_4$Autres[9:12], Tab3_4$`Erreur de Monte Carlo`[12],
  Tab3_4$Autres[13:16], Tab3_4$`Erreur de Monte Carlo`[16]
  
)


# nsim = 1000, Scénarios 01, 02, 03 et 04

comparaison01$`Scénario01_0%` <- sapply(comparaison01$`Scénario01_0%`, FUN = function(x)x[[2]])
comparaison01$`Scénario02_10%` <- sapply(comparaison01$`Scénario02_10%`, FUN = function(x)x[[2]])
comparaison01$`Scénario03_20%` <- sapply(comparaison01$`Scénario03_20%`, FUN = function(x)x[[2]])
comparaison01$`Scénario04_40%` <- sapply(comparaison01$`Scénario04_40%`, FUN = function(x)x[[2]])

View(comparaison01)
################################################################################
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario05_10%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario05_15%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario05_30%.RData")

# nsim = 1000, Scénarios 05 prévalence de I2  ~ 10%, 15% et 30%, taux de co-infection ~30%-40%

comparaison02 <- data.frame(matrix(ncol = 5, 
                                   nrow = 30))

colnames(comparaison02) <- c(
  "Estimation", 
  "Performance",
  "Scénario05_10%",
  "Scénario05_15%",
  "Scénario05_30%"
)

comparaison02$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison02$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison02$`Scénario05_10%` <- c(
    
    Tab_5_1$Autres, Tab_5_1$`Erreur de Monte Carlo`[4],
    Tab2_5_1$Autres, Tab2_5_1$`Erreur de Monte Carlo`[4],
    Tab3_5_1$Autres[1:4], Tab3_5_1$`Erreur de Monte Carlo`[4],
    Tab3_5_1$Autres[5:8], Tab3_5_1$`Erreur de Monte Carlo`[8],
    Tab3_5_1$Autres[9:12], Tab3_5_1$`Erreur de Monte Carlo`[12],
    Tab3_5_1$Autres[13:16], Tab3_5_1$`Erreur de Monte Carlo`[16]
  
)

comparaison02$`Scénario05_15%` <- c(
  
  Tab_5_2$Autres, Tab_5_2$`Erreur de Monte Carlo`[4],
  Tab2_5_2$Autres, Tab2_5_2$`Erreur de Monte Carlo`[4],
  Tab3_5_2$Autres[1:4], Tab3_5_2$`Erreur de Monte Carlo`[4],
  Tab3_5_2$Autres[5:8], Tab3_5_2$`Erreur de Monte Carlo`[8],
  Tab3_5_2$Autres[9:12], Tab3_5_2$`Erreur de Monte Carlo`[12],
  Tab3_5_2$Autres[13:16], Tab3_5_2$`Erreur de Monte Carlo`[16]
  
)

comparaison02$`Scénario05_30%` <- c(
  
  Tab_5_3$Autres, Tab_5_3$`Erreur de Monte Carlo`[4],
  Tab2_5_3$Autres, Tab2_5_3$`Erreur de Monte Carlo`[4],
  Tab3_5_3$Autres[1:4], Tab3_5_3$`Erreur de Monte Carlo`[4],
  Tab3_5_3$Autres[5:8], Tab3_5_3$`Erreur de Monte Carlo`[8],
  Tab3_5_3$Autres[9:12], Tab3_5_3$`Erreur de Monte Carlo`[12],
  Tab3_5_3$Autres[13:16], Tab3_5_3$`Erreur de Monte Carlo`[16]
  
)

comparaison02$`Scénario05_10%` <- sapply(comparaison02$`Scénario05_10%`, FUN = function(x)x[[2]])
comparaison02$`Scénario05_15%` <- sapply(comparaison02$`Scénario05_15%`, FUN = function(x)x[[2]])
comparaison02$`Scénario05_30%` <- sapply(comparaison02$`Scénario05_30%`, FUN = function(x)x[[2]])

View(comparaison02)
################################################################################
####################   Représentations graphiques  #############################
################################################################################
## Scénarios avec différents taux de co-infection

dat <- data.frame(x_biais = c(comparaison01$`Scénario01_0%`[[1]], comparaison01$`Scénario02_10%`[[1]], comparaison01$`Scénario03_20%`[[1]], comparaison01$`Scénario04_40%`[[1]],
                   comparaison01$`Scénario01_0%`[[6]], comparaison01$`Scénario02_10%`[[6]], comparaison01$`Scénario03_20%`[[6]], comparaison01$`Scénario04_40%`[[6]],
                   comparaison01$`Scénario01_0%`[[11]], comparaison01$`Scénario02_10%`[[11]], comparaison01$`Scénario03_20%`[[11]], comparaison01$`Scénario04_40%`[[11]],
                  comparaison01$`Scénario01_0%`[[16]], comparaison01$`Scénario02_10%`[[16]], comparaison01$`Scénario03_20%`[[16]], comparaison01$`Scénario04_40%`[[16]],
                  comparaison01$`Scénario01_0%`[[21]], comparaison01$`Scénario02_10%`[[21]], comparaison01$`Scénario03_20%`[[21]], comparaison01$`Scénario04_40%`[[21]],
                  comparaison01$`Scénario01_0%`[[26]], comparaison01$`Scénario02_10%`[[26]], comparaison01$`Scénario03_20%`[[26]], comparaison01$`Scénario04_40%`[[26]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 4),
                  Scénario = c(1, 2, 3, 4))
                               

dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)

 
ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

## Scénarios avec différentes valeurs de prévalence de I2 et un taux de co-infection entre 30% et 40%

dat1 <- data.frame(x_biais = c(comparaison02$`Scénario05_10%`[[1]], comparaison02$`Scénario05_15%`[[1]], comparaison02$`Scénario05_30%`[[1]],
                              comparaison02$`Scénario05_10%`[[6]], comparaison02$`Scénario05_15%`[[6]], comparaison02$`Scénario05_30%`[[6]],
                              comparaison02$`Scénario05_10%`[[11]], comparaison02$`Scénario05_15%`[[11]], comparaison02$`Scénario05_30%`[[11]],
                              comparaison02$`Scénario05_10%`[[16]], comparaison02$`Scénario05_15%`[[16]], comparaison02$`Scénario05_30%`[[16]],
                              comparaison02$`Scénario05_10%`[[21]], comparaison02$`Scénario05_15%`[[21]], comparaison02$`Scénario05_30%`[[21]],
                              comparaison02$`Scénario05_10%`[[26]], comparaison02$`Scénario05_15%`[[26]], comparaison02$`Scénario05_30%`[[26]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                  Scénario = c(1, 2, 3))


dat1$Scénario <- as.factor(dat1$Scénario)
dat1$Méthode<- as.factor(dat1$Méthode)


ggplot(dat1) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))
