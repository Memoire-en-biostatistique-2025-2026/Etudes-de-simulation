# Chargement des librairies nécessaires

library(ggplot2)
library(patchwork)
library(kableExtra)
library("gridExtra")
library("cowplot")
library(ggpubr)
library(Rmisc)

################################################################################
## Taille d'échantillo = 1000
################################################################################
## Scénarios avec différents taux de co-infection
# Chargement des résultats bruts

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
####################### Représentation graphique ###############################

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
################################################################################
# Scénarios avec différentes valeurs de prévalence de I2 et un taux de co-infection entre 30% et 40%

## Chargement des résultats bruts

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

####################### Représentation graphique ###############################

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
################################################################################
# Scénarios avec différentes valeurs de couverture vaccinale et un taux de co-infection entre 30% et 40%

## Chargement des résultats bruts

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario06_15%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario06_50%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario06_70%.RData")

# nsim = 1000, Scénarios 06 couverture vaccinale  ~ 15%, 50% et 70%, taux de co-infection ~30%-40%

comparaison03 <- data.frame(matrix(ncol = 5, 
                                   nrow = 30))

colnames(comparaison03) <- c(
  "Estimation", 
  "Performance",
  "Scénario06_15%",
  "Scénario06_50%",
  "Scénario06_70%"
)

comparaison03$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison03$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison03$`Scénario06_15%` <- c(
  
  Tab_6_1$Autres, Tab_6_1$`Erreur de Monte Carlo`[4],
  Tab2_6_1$Autres, Tab2_6_1$`Erreur de Monte Carlo`[4],
  Tab3_6_1$Autres[1:4], Tab3_6_1$`Erreur de Monte Carlo`[4],
  Tab3_6_1$Autres[5:8], Tab3_6_1$`Erreur de Monte Carlo`[8],
  Tab3_6_1$Autres[9:12], Tab3_6_1$`Erreur de Monte Carlo`[12],
  Tab3_6_1$Autres[13:16], Tab3_6_1$`Erreur de Monte Carlo`[16]
  
)

comparaison03$`Scénario06_50%` <- c(
  
  Tab_6_2$Autres, Tab_6_2$`Erreur de Monte Carlo`[4],
  Tab2_6_2$Autres, Tab2_6_2$`Erreur de Monte Carlo`[4],
  Tab3_6_2$Autres[1:4], Tab3_6_2$`Erreur de Monte Carlo`[4],
  Tab3_6_2$Autres[5:8], Tab3_6_2$`Erreur de Monte Carlo`[8],
  Tab3_6_2$Autres[9:12], Tab3_6_2$`Erreur de Monte Carlo`[12],
  Tab3_6_2$Autres[13:16], Tab3_6_2$`Erreur de Monte Carlo`[16]
  
)

comparaison03$`Scénario06_70%` <- c(
  
  Tab_6_3$Autres, Tab_6_3$`Erreur de Monte Carlo`[4],
  Tab2_6_3$Autres, Tab2_6_3$`Erreur de Monte Carlo`[4],
  Tab3_6_3$Autres[1:4], Tab3_6_3$`Erreur de Monte Carlo`[4],
  Tab3_6_3$Autres[5:8], Tab3_6_3$`Erreur de Monte Carlo`[8],
  Tab3_6_3$Autres[9:12], Tab3_6_3$`Erreur de Monte Carlo`[12],
  Tab3_6_3$Autres[13:16], Tab3_6_3$`Erreur de Monte Carlo`[16]
  
)

comparaison03$`Scénario06_15%` <- sapply(comparaison03$`Scénario06_15%`, FUN = function(x)x[[2]])
comparaison03$`Scénario06_50%` <- sapply(comparaison03$`Scénario06_50%`, FUN = function(x)x[[2]])
comparaison03$`Scénario06_70%` <- sapply(comparaison03$`Scénario06_70%`, FUN = function(x)x[[2]])

View(comparaison03)

####################### Représentation graphique ###############################

dat2 <- data.frame(x_biais = c(comparaison03$`Scénario06_15%`[[1]], comparaison03$`Scénario06_50%`[[1]], comparaison03$`Scénario06_70%`[[1]],
                               comparaison03$`Scénario06_15%`[[6]], comparaison03$`Scénario06_50%`[[6]], comparaison03$`Scénario06_70%`[[6]],
                               comparaison03$`Scénario06_15%`[[11]], comparaison03$`Scénario06_50%`[[11]], comparaison03$`Scénario06_70%`[[11]],
                               comparaison03$`Scénario06_15%`[[16]], comparaison03$`Scénario06_50%`[[16]], comparaison03$`Scénario06_70%`[[16]],
                               comparaison03$`Scénario06_15%`[[21]], comparaison03$`Scénario06_50%`[[21]], comparaison03$`Scénario06_70%`[[21]],
                               comparaison03$`Scénario06_15%`[[26]], comparaison03$`Scénario06_50%`[[26]], comparaison03$`Scénario06_70%`[[26]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                   Scénario = c(1, 2, 3))


dat2$Scénario <- as.factor(dat2$Scénario)
dat2$Méthode<- as.factor(dat2$Méthode)


ggplot(dat2) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

################################################################################
# Scénarios avec différentes valeurs de prévalence de I1 et un taux de co-infection entre 30% et 40%

## Chargement des résultats bruts

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario07_10%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario07_30%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario07_50%.RData")

# nsim = 1000, Scénarios 07 prévalence de I1  ~ 10%, 30% et 50%, taux de co-infection ~30%-40%

comparaison04 <- data.frame(matrix(ncol = 5, 
                                   nrow = 30))

colnames(comparaison04) <- c(
  "Estimation", 
  "Performance",
  "Scénario07_10%",
  "Scénario07_30%",
  "Scénario07_50%"
)

comparaison04$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison04$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison04$`Scénario07_10%` <- c(
  
  Tab_7_1$Autres, Tab_7_1$`Erreur de Monte Carlo`[4],
  Tab2_7_1$Autres, Tab2_7_1$`Erreur de Monte Carlo`[4],
  Tab3_7_1$Autres[1:4], Tab3_7_1$`Erreur de Monte Carlo`[4],
  Tab3_7_1$Autres[5:8], Tab3_7_1$`Erreur de Monte Carlo`[8],
  Tab3_7_1$Autres[9:12], Tab3_7_1$`Erreur de Monte Carlo`[12],
  Tab3_7_1$Autres[13:16], Tab3_7_1$`Erreur de Monte Carlo`[16]
  
)

comparaison04$`Scénario07_30%` <- c(
  
  Tab_7_2$Autres, Tab_7_2$`Erreur de Monte Carlo`[4],
  Tab2_7_2$Autres, Tab2_7_2$`Erreur de Monte Carlo`[4],
  Tab3_7_2$Autres[1:4], Tab3_7_2$`Erreur de Monte Carlo`[4],
  Tab3_7_2$Autres[5:8], Tab3_7_2$`Erreur de Monte Carlo`[8],
  Tab3_7_2$Autres[9:12], Tab3_7_2$`Erreur de Monte Carlo`[12],
  Tab3_7_2$Autres[13:16], Tab3_7_2$`Erreur de Monte Carlo`[16]
  
)

comparaison04$`Scénario07_50%` <- c(
  
  Tab_7_3$Autres, Tab_7_3$`Erreur de Monte Carlo`[4],
  Tab2_7_3$Autres, Tab2_7_3$`Erreur de Monte Carlo`[4],
  Tab3_7_3$Autres[1:4], Tab3_7_3$`Erreur de Monte Carlo`[4],
  Tab3_7_3$Autres[5:8], Tab3_7_3$`Erreur de Monte Carlo`[8],
  Tab3_7_3$Autres[9:12], Tab3_7_3$`Erreur de Monte Carlo`[12],
  Tab3_7_3$Autres[13:16], Tab3_7_3$`Erreur de Monte Carlo`[16]
  
)

comparaison04$`Scénario07_10%` <- sapply(comparaison04$`Scénario07_10%`, FUN = function(x)x[[2]])
comparaison04$`Scénario07_30%` <- sapply(comparaison04$`Scénario07_30%`, FUN = function(x)x[[2]])
comparaison04$`Scénario07_50%` <- sapply(comparaison04$`Scénario07_50%`, FUN = function(x)x[[2]])

View(comparaison04)

####################### Représentation graphique ###############################

dat3 <- data.frame(x_biais = c(comparaison04$`Scénario07_10%`[[1]], comparaison04$`Scénario07_30%`[[1]], comparaison04$`Scénario07_50%`[[1]],
                               comparaison04$`Scénario07_10%`[[6]], comparaison04$`Scénario07_30%`[[6]], comparaison04$`Scénario07_50%`[[6]],
                               comparaison04$`Scénario07_10%`[[11]], comparaison04$`Scénario07_30%`[[11]], comparaison04$`Scénario07_50%`[[11]],
                               comparaison04$`Scénario07_10%`[[16]], comparaison04$`Scénario07_30%`[[16]], comparaison04$`Scénario07_50%`[[16]],
                               comparaison04$`Scénario07_10%`[[21]], comparaison04$`Scénario07_30%`[[21]], comparaison04$`Scénario07_50%`[[21]],
                               comparaison04$`Scénario07_10%`[[26]], comparaison04$`Scénario07_30%`[[26]], comparaison04$`Scénario07_50%`[[26]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                   Scénario = c(1, 2, 3))


dat3$Scénario <- as.factor(dat3$Scénario)
dat3$Méthode<- as.factor(dat3$Méthode)


ggplot(dat3) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

################################################################################
## Taille d'échantillo = 5000
################################################################################

# Tableau de comparaison 

comparaison05 <- data.frame(matrix(ncol = 6, 
                                   nrow = 30))

colnames(comparaison05) <- c(
  "Estimation", 
  "Performance",
  "Scénario01_0%",
  "Scénario02_10%",
  "Scénario03_20%",
  "Scénario04_40%"
)

comparaison05$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison05$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison05$`Scénario01_0%` <- c(
  
  -0.000521146215135881, 0.00141551900341267, 0.00141437507778681, 0.0192973888418796, 0.877,
  -0.00261270050723522, 0.00149261895418944, 0.00141437507778681, 0.0234009034199065, 0.865,
  -0.00281469760172037,  0.00162864036378977, 0.00163493424601511, 0.0236851250420008, 0.856,
  0.002428957699253, 0.00157811936632933, 0.00158244108246776, 0.049307571941646, 0.634,
  -0.00780576867498711, 0.00476654458999823, 0.00482270807001564, 0.377731629823567, 0.884,
  -0.00840737636418898, 0.00294625647934802, 0.00301399420019779, 0.385298359516229, 0.873
  
)

comparaison05$`Scénario02_10%` <- c(
  
  0.0346021559456765, 0.00169746059235832, 0.00289307232785488, 0.0202695772601545, 0.752,   
  0.0266241900955686, 0.00178820138074282, 0.00289307232785488, 0.0231930872417901, 0.793,
  0.0262127431787281, 0.0019700580050842, 0.00265519585203307, 0.0253505303389045, 0.831,
  0.0321714432263664, 0.00190814254587139, 0.00294123616259284, 0.0560217757419352 , 0.49,
  0.0153256782006115, 0.00720737672879046, 0.00743504576437037, 0.597017343022217, 0.874,
  0.0124125795822325, 0.00755634402191734, 0.00770285980978068,  0.794091658462609, 0.86
  
)

comparaison05$`Scénario03_20%` <- c(
  
   0.0702207334602339, 0.00203205380131056, 0.00696097315520247, 0.018004067803242, 0.483,
   0.0547388990180326, 0.00203577355422482, 0.00696097315520247, 0.0169013017245136, 0.637,
   0.0527365698189655, 0.00227903497533604, 0.00505790173663133, 0.0230005410149488, 0.685,
   0.0619079688725061, 0.00266711719755541, 0.00649704669027704, 0.0555927403366925, 0.361,
   0.0376474201503548, 0.00925173035744706, 0.010659806871067, 0.615739002376905, 0.798,
   0.0285351902425544, 0.013939257389593, 0.0147395752143822, 1.12084275720326, 0.772 
  
)

comparaison05$`Scénario04_40%` <- c(
  
  0.144361215601368, 0.00640376926685333 , 0.0272375260674911, 0.0291408824851266, 0.319,
  0.119883432905491, 0.00661862847792735, 0.0272375260674911, 0.0224086167200437, 0.528,
  0.114389824802602, 0.00741330070480691, 0.0204909194224722, 0.0365185464025841, 0.638,
  -0.00473479892132594, 0.0617660214527524, 0.0617266737521251, 1.91949706896884, 0.99,
  0.077570090978503, 0.0310872524225409, 0.0370732841845316, 0.920424519359417, 0.794,
  -0.00146456469367395, 0.0678527259288813, 0.0677870181526944, 1.43411748434497, 0.803
  
)

View(comparaison05)

dat4 <- data.frame(x_biais = c(comparaison05$`Scénario01_0%`[[1]], comparaison05$`Scénario02_10%`[[1]], comparaison05$`Scénario03_20%`[[1]], comparaison05$`Scénario04_40%`[[1]],
                              comparaison05$`Scénario01_0%`[[6]], comparaison05$`Scénario02_10%`[[6]], comparaison05$`Scénario03_20%`[[6]], comparaison05$`Scénario04_40%`[[6]],
                              comparaison05$`Scénario01_0%`[[11]], comparaison05$`Scénario02_10%`[[11]], comparaison05$`Scénario03_20%`[[11]], comparaison05$`Scénario04_40%`[[11]],
                              comparaison05$`Scénario01_0%`[[16]], comparaison05$`Scénario02_10%`[[16]], comparaison05$`Scénario03_20%`[[16]], comparaison05$`Scénario04_40%`[[16]],
                              comparaison05$`Scénario01_0%`[[21]], comparaison05$`Scénario02_10%`[[21]], comparaison05$`Scénario03_20%`[[21]], comparaison05$`Scénario04_40%`[[21]],
                              comparaison05$`Scénario01_0%`[[26]], comparaison05$`Scénario02_10%`[[26]], comparaison05$`Scénario03_20%`[[26]], comparaison05$`Scénario04_40%`[[26]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 4),
                  Scénario = c(1, 2, 3, 4))


dat4$Scénario <- as.factor(dat4$Scénario)
dat4$Méthode<- as.factor(dat4$Méthode)


ggplot(dat4) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))
