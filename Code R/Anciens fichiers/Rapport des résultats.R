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
## Scénarios avec différents taux de co-infection : 0%, 10%, 20% et 40%
# nsim = 1000, Scénarios 01, 02, 03 et 04

# Chargement des résultats bruts

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario01.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario02.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario03.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario04.RData")

# Tableau de comparaison 

comparaison01 <- data.frame(matrix(ncol = 6, 
                                   nrow = 42))

colnames(comparaison01) <- c(
  "Estimation", 
  "Performance",
  "Scénario01_0%",
  "Scénario02_10%",
  "Scénario03_20%",
  "Scénario04_40%"
)

comparaison01$Estimation <- c("RegLog", "-", "-", "-", "-","-", "IPW", "-", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-", "-",
                              "TNDDR_GLM", "-", "-", "-", "-", "-")

comparaison01$Performance <- rep(c("Biais_med", "Biais_moy", "Variance", "MSE", "Précision", "%Cov"), 7)

comparaison01$`Scénario01_0%` <- c(
  
  (median(resultats_1$RRc) - mean(l_vraiRRc_1)), Tab01_1$Autres, Tab01_1$`Erreur de Monte Carlo`[4],
  (median(resultats2_1$RRm) - mean(l_vraiRRm_1)), Tab02_1$Autres, Tab02_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_1$RRm_RF) - mean(l_vraiRRm_1)), Tab03_1$Autres[1:4], Tab03_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_1$RRm_Lasso) - mean(l_vraiRRm_1)),Tab03_1$Autres[5:8], Tab03_1$`Erreur de Monte Carlo`[8],
  (median(resultats3_1$RRm_Mars) - mean(l_vraiRRm_1)), Tab03_1$Autres[9:12], Tab03_1$`Erreur de Monte Carlo`[12],
  (median(resultats3_1$RRm_RN) - mean(l_vraiRRm_1)), Tab03_1$Autres[13:16], Tab03_1$`Erreur de Monte Carlo`[16],
  (median(resultats3_1$RRm_GLM) - mean(l_vraiRRm_1)), Tab03_1$Autres[17:20], Tab03_1$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario02_10%` <- c(
  
  (median(resultats_2$RRc) - mean(l_vraiRRc_2)), Tab01_2$Autres, Tab01_2$`Erreur de Monte Carlo`[4],
  (median(resultats2_2$RRm) - mean(l_vraiRRm_2)), Tab02_2$Autres, Tab02_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_2$RRm_RF) - mean(l_vraiRRm_2)), Tab03_2$Autres[1:4], Tab03_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_2$RRm_Lasso) - mean(l_vraiRRm_2)),Tab03_2$Autres[5:8], Tab03_2$`Erreur de Monte Carlo`[8],
  (median(resultats3_2$RRm_Mars) - mean(l_vraiRRm_2)), Tab03_2$Autres[9:12], Tab03_2$`Erreur de Monte Carlo`[12],
  (median(resultats3_2$RRm_RN) - mean(l_vraiRRm_2)), Tab03_2$Autres[13:16], Tab03_2$`Erreur de Monte Carlo`[16],
  (median(resultats3_2$RRm_GLM) - mean(l_vraiRRm_2)), Tab03_2$Autres[17:20], Tab03_2$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario03_20%` <- c(
  
  (median(resultats_3$RRc) - mean(l_vraiRRc_3)), Tab01_3$Autres, Tab01_3$`Erreur de Monte Carlo`[4],
  (median(resultats2_3$RRm) - mean(l_vraiRRm_3)), Tab02_3$Autres, Tab02_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_3$RRm_RF) - mean(l_vraiRRm_3)), Tab03_3$Autres[1:4], Tab03_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_3$RRm_Lasso) - mean(l_vraiRRm_3)),Tab03_3$Autres[5:8], Tab03_3$`Erreur de Monte Carlo`[8],
  (median(resultats3_3$RRm_Mars) - mean(l_vraiRRm_3)), Tab03_3$Autres[9:12], Tab03_3$`Erreur de Monte Carlo`[12],
  (median(resultats3_3$RRm_RN) - mean(l_vraiRRm_3)), Tab03_3$Autres[13:16], Tab03_3$`Erreur de Monte Carlo`[16],
  (median(resultats3_3$RRm_GLM) - mean(l_vraiRRm_3)), Tab03_3$Autres[17:20], Tab03_3$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario04_40%` <- c(
  
  (median(resultats_4$RRc) - mean(l_vraiRRc_4)), Tab01_4$Autres, Tab01_4$`Erreur de Monte Carlo`[4],
  (median(resultats2_4$RRm) - mean(l_vraiRRm_4)), Tab02_4$Autres, Tab02_4$`Erreur de Monte Carlo`[4],
  (median(resultats3_4$RRm_RF) - mean(l_vraiRRm_4)), Tab03_4$Autres[1:4], Tab03_4$`Erreur de Monte Carlo`[4],
  (median(resultats3_4$RRm_Lasso, na.rm = TRUE) - mean(l_vraiRRm_4)),Tab03_4$Autres[5:8], Tab03_4$`Erreur de Monte Carlo`[8],
  (median(resultats3_4$RRm_Mars) - mean(l_vraiRRm_4)), Tab03_4$Autres[9:12], Tab03_4$`Erreur de Monte Carlo`[12],
  (median(resultats3_4$RRm_RN) - mean(l_vraiRRm_4)), Tab03_4$Autres[13:16], Tab03_4$`Erreur de Monte Carlo`[16],
  (median(resultats3_4$RRm_GLM) - mean(l_vraiRRm_4)), Tab03_4$Autres[17:20], Tab03_4$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario01_0%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario01_0%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison01$`Scénario02_10%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario02_10%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison01$`Scénario03_20%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario03_20%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison01$`Scénario04_40%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario04_40%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison01)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison01$`Scénario01_0%`[[2]], comparaison01$`Scénario02_10%`[[2]], comparaison01$`Scénario03_20%`[[2]], comparaison01$`Scénario04_40%`[[2]],
                              comparaison01$`Scénario01_0%`[[8]], comparaison01$`Scénario02_10%`[[8]], comparaison01$`Scénario03_20%`[[8]], comparaison01$`Scénario04_40%`[[8]],
                              comparaison01$`Scénario01_0%`[[14]], comparaison01$`Scénario02_10%`[[14]], comparaison01$`Scénario03_20%`[[14]], comparaison01$`Scénario04_40%`[[14]],
                              comparaison01$`Scénario01_0%`[[20]], comparaison01$`Scénario02_10%`[[20]], comparaison01$`Scénario03_20%`[[20]], comparaison01$`Scénario04_40%`[[20]],
                              comparaison01$`Scénario01_0%`[[26]], comparaison01$`Scénario02_10%`[[26]], comparaison01$`Scénario03_20%`[[26]], comparaison01$`Scénario04_40%`[[26]],
                              comparaison01$`Scénario01_0%`[[32]], comparaison01$`Scénario02_10%`[[32]], comparaison01$`Scénario03_20%`[[32]], comparaison01$`Scénario04_40%`[[32]],
                              comparaison01$`Scénario01_0%`[[38]], comparaison01$`Scénario02_10%`[[38]], comparaison01$`Scénario03_20%`[[38]], comparaison01$`Scénario04_40%`[[38]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 4),
                  Scénario = c(1, 2, 3, 4))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture

dat <- data.frame(x_cov = c(comparaison01$`Scénario01_0%`[[6]], comparaison01$`Scénario02_10%`[[6]], comparaison01$`Scénario03_20%`[[6]], comparaison01$`Scénario04_40%`[[6]],
                              comparaison01$`Scénario01_0%`[[12]], comparaison01$`Scénario02_10%`[[12]], comparaison01$`Scénario03_20%`[[12]], comparaison01$`Scénario04_40%`[[12]],
                              comparaison01$`Scénario01_0%`[[18]], comparaison01$`Scénario02_10%`[[18]], comparaison01$`Scénario03_20%`[[18]], comparaison01$`Scénario04_40%`[[18]],
                              comparaison01$`Scénario01_0%`[[24]], comparaison01$`Scénario02_10%`[[24]], comparaison01$`Scénario03_20%`[[24]], comparaison01$`Scénario04_40%`[[24]],
                              comparaison01$`Scénario01_0%`[[30]], comparaison01$`Scénario02_10%`[[30]], comparaison01$`Scénario03_20%`[[30]], comparaison01$`Scénario04_40%`[[30]],
                              comparaison01$`Scénario01_0%`[[36]], comparaison01$`Scénario02_10%`[[36]], comparaison01$`Scénario03_20%`[[36]], comparaison01$`Scénario04_40%`[[36]],
                              comparaison01$`Scénario01_0%`[[42]], comparaison01$`Scénario02_10%`[[42]], comparaison01$`Scénario03_20%`[[42]], comparaison01$`Scénario04_40%`[[42]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 4),
                  Scénario = c(1, 2, 3, 4))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Biais médian

dat <- data.frame(x_biais = c(comparaison01$`Scénario01_0%`[[1]], comparaison01$`Scénario02_10%`[[1]], comparaison01$`Scénario03_20%`[[1]], comparaison01$`Scénario04_40%`[[1]],
                              comparaison01$`Scénario01_0%`[[7]], comparaison01$`Scénario02_10%`[[7]], comparaison01$`Scénario03_20%`[[7]], comparaison01$`Scénario04_40%`[[7]],
                              comparaison01$`Scénario01_0%`[[13]], comparaison01$`Scénario02_10%`[[13]], comparaison01$`Scénario03_20%`[[13]], comparaison01$`Scénario04_40%`[[13]],
                              comparaison01$`Scénario01_0%`[[19]], comparaison01$`Scénario02_10%`[[19]], comparaison01$`Scénario03_20%`[[19]], comparaison01$`Scénario04_40%`[[19]],
                              comparaison01$`Scénario01_0%`[[25]], comparaison01$`Scénario02_10%`[[25]], comparaison01$`Scénario03_20%`[[25]], comparaison01$`Scénario04_40%`[[25]],
                              comparaison01$`Scénario01_0%`[[31]], comparaison01$`Scénario02_10%`[[31]], comparaison01$`Scénario03_20%`[[31]], comparaison01$`Scénario04_40%`[[31]],
                              comparaison01$`Scénario01_0%`[[37]], comparaison01$`Scénario02_10%`[[37]], comparaison01$`Scénario03_20%`[[37]], comparaison01$`Scénario04_40%`[[37]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 4),
                  Scénario = c(1, 2, 3, 4))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

################################################################################
################################################################################
# Scénarios avec différentes valeurs de prévalence de I2 et un taux de co-infection entre 30% et 40%

## Chargement des résultats bruts

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario05_10%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario05_15%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario05_30%.RData")

# nsim = 1000, Scénarios 05 prévalence de I2  ~ 10%, 15% et 30%, taux de co-infection ~30%-40%

comparaison02 <- data.frame(matrix(ncol = 5, 
                                   nrow = 42))

colnames(comparaison02) <- c(
  "Estimation", 
  "Performance",
  "Scénario05_10%",
  "Scénario05_15%",
  "Scénario05_30%"
)

comparaison02$Estimation <- c("RegLog", "-", "-", "-", "-","-", "IPW", "-", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-", "-",
                              "TNDDR_GLM", "-", "-", "-", "-", "-")

comparaison02$Performance <- rep(c("Biais_med", "Biais_moy", "Variance", "MSE", "Précision", "%Cov"), 7)

comparaison02$`Scénario05_10%` <- c(
  
  (median(resultats_5_1$RRc) - mean(l_vraiRRc_5_1)), Tab01_5_1$Autres, Tab01_5_1$`Erreur de Monte Carlo`[4],
  (median(resultats2_5_1$RRm) - mean(l_vraiRRm_5_1)), Tab02_5_1$Autres, Tab02_5_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_1$RRm_RF) - mean(l_vraiRRm_5_1)), Tab03_5_1$Autres[1:4], Tab03_5_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_1$RRm_Lasso) - mean(l_vraiRRm_5_1)),Tab03_5_1$Autres[5:8], Tab03_5_1$`Erreur de Monte Carlo`[8],
  (median(resultats3_5_1$RRm_Mars) - mean(l_vraiRRm_5_1)), Tab03_5_1$Autres[9:12], Tab03_5_1$`Erreur de Monte Carlo`[12],
  (median(resultats3_5_1$RRm_RN) - mean(l_vraiRRm_5_1)), Tab03_5_1$Autres[13:16], Tab03_5_1$`Erreur de Monte Carlo`[16],
  (median(resultats3_5_1$RRm_GLM) - mean(l_vraiRRm_5_1)), Tab03_5_1$Autres[17:20], Tab03_5_1$`Erreur de Monte Carlo`[20]
  
)

comparaison02$`Scénario05_15%` <- c(
  
  (median(resultats_5_2$RRc) - mean(l_vraiRRc_5_2)), Tab01_5_2$Autres, Tab01_5_2$`Erreur de Monte Carlo`[4],
  (median(resultats2_5_2$RRm) - mean(l_vraiRRm_5_2)), Tab02_5_2$Autres, Tab02_5_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_2$RRm_RF) - mean(l_vraiRRm_5_2)), Tab03_5_2$Autres[1:4], Tab03_5_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_2$RRm_Lasso) - mean(l_vraiRRm_5_2)),Tab03_5_2$Autres[5:8], Tab03_5_2$`Erreur de Monte Carlo`[8],
  (median(resultats3_5_2$RRm_Mars) - mean(l_vraiRRm_5_2)), Tab03_5_2$Autres[9:12], Tab03_5_2$`Erreur de Monte Carlo`[12],
  (median(resultats3_5_2$RRm_RN) - mean(l_vraiRRm_5_2)), Tab03_5_2$Autres[13:16], Tab03_5_2$`Erreur de Monte Carlo`[16],
  (median(resultats3_5_2$RRm_GLM) - mean(l_vraiRRm_5_2)), Tab03_5_2$Autres[17:20], Tab03_5_2$`Erreur de Monte Carlo`[20]
  
)

comparaison02$`Scénario05_30%` <- c(
  
  (median(resultats_5_3$RRc) - mean(l_vraiRRc_5_3)), Tab01_5_3$Autres, Tab01_5_3$`Erreur de Monte Carlo`[4],
  (median(resultats2_5_3$RRm) - mean(l_vraiRRm_5_3)), Tab02_5_3$Autres, Tab02_5_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_3$RRm_RF) - mean(l_vraiRRm_5_3)), Tab03_5_3$Autres[1:4], Tab03_5_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_3$RRm_Lasso) - mean(l_vraiRRm_5_3)),Tab03_5_3$Autres[5:8], Tab03_5_3$`Erreur de Monte Carlo`[8],
  (median(resultats3_5_3$RRm_Mars) - mean(l_vraiRRm_5_3)), Tab03_5_3$Autres[9:12], Tab03_5_3$`Erreur de Monte Carlo`[12],
  (median(resultats3_5_3$RRm_RN) - mean(l_vraiRRm_5_3)), Tab03_5_3$Autres[13:16], Tab03_5_3$`Erreur de Monte Carlo`[16],
  (median(resultats3_5_3$RRm_GLM) - mean(l_vraiRRm_5_3)), Tab03_5_3$Autres[17:20], Tab03_5_3$`Erreur de Monte Carlo`[20]
  
)

comparaison02$`Scénario05_10%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison02$`Scénario05_10%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison02$`Scénario05_15%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison02$`Scénario05_15%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison02$`Scénario05_30%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison02$`Scénario05_30%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison02)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison02$`Scénario05_10%`[[2]], comparaison02$`Scénario05_15%`[[2]], comparaison02$`Scénario05_30%`[[2]],
                               comparaison02$`Scénario05_10%`[[8]], comparaison02$`Scénario05_15%`[[8]], comparaison02$`Scénario05_30%`[[8]],
                               comparaison02$`Scénario05_10%`[[14]], comparaison02$`Scénario05_15%`[[14]], comparaison02$`Scénario05_30%`[[14]],
                               comparaison02$`Scénario05_10%`[[20]], comparaison02$`Scénario05_15%`[[20]], comparaison02$`Scénario05_30%`[[20]],
                               comparaison02$`Scénario05_10%`[[26]], comparaison02$`Scénario05_15%`[[26]], comparaison02$`Scénario05_30%`[[26]],
                               comparaison02$`Scénario05_10%`[[32]], comparaison02$`Scénario05_15%`[[32]], comparaison02$`Scénario05_30%`[[32]],
                              comparaison02$`Scénario05_10%`[[38]], comparaison02$`Scénario05_15%`[[38]], comparaison02$`Scénario05_30%`[[38]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                   Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture

dat <- data.frame(x_cov = c(comparaison02$`Scénario05_10%`[[6]], comparaison02$`Scénario05_15%`[[6]], comparaison02$`Scénario05_30%`[[6]],
                               comparaison02$`Scénario05_10%`[[12]], comparaison02$`Scénario05_15%`[[12]], comparaison02$`Scénario05_30%`[[12]],
                               comparaison02$`Scénario05_10%`[[18]], comparaison02$`Scénario05_15%`[[18]], comparaison02$`Scénario05_30%`[[18]],
                               comparaison02$`Scénario05_10%`[[24]], comparaison02$`Scénario05_15%`[[24]], comparaison02$`Scénario05_30%`[[24]],
                               comparaison02$`Scénario05_10%`[[30]], comparaison02$`Scénario05_15%`[[30]], comparaison02$`Scénario05_30%`[[30]],
                               comparaison02$`Scénario05_10%`[[36]], comparaison02$`Scénario05_15%`[[36]], comparaison02$`Scénario05_30%`[[36]],
                               comparaison02$`Scénario05_10%`[[42]], comparaison02$`Scénario05_15%`[[42]], comparaison02$`Scénario05_30%`[[42]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                   Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Biais median

dat <- data.frame(x_biais = c(comparaison02$`Scénario05_10%`[[1]], comparaison02$`Scénario05_15%`[[1]], comparaison02$`Scénario05_30%`[[1]],
                              comparaison02$`Scénario05_10%`[[7]], comparaison02$`Scénario05_15%`[[7]], comparaison02$`Scénario05_30%`[[7]],
                              comparaison02$`Scénario05_10%`[[13]], comparaison02$`Scénario05_15%`[[13]], comparaison02$`Scénario05_30%`[[13]],
                              comparaison02$`Scénario05_10%`[[19]], comparaison02$`Scénario05_15%`[[19]], comparaison02$`Scénario05_30%`[[19]],
                              comparaison02$`Scénario05_10%`[[25]], comparaison02$`Scénario05_15%`[[25]], comparaison02$`Scénario05_30%`[[25]],
                              comparaison02$`Scénario05_10%`[[31]], comparaison02$`Scénario05_15%`[[31]], comparaison02$`Scénario05_30%`[[31]],
                              comparaison02$`Scénario05_10%`[[37]], comparaison02$`Scénario05_15%`[[37]], comparaison02$`Scénario05_30%`[[37]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

################################################################################
################################################################################
# Scénarios avec différentes valeurs de couverture vaccinale et un taux de co-infection entre 30% et 40%

## Chargement des résultats bruts

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario06_50%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario06_70%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario06_85%.RData")

# nsim = 1000, Scénarios 06 couverture vaccinale  ~ 50%, 70% et 85%, taux de co-infection ~30%-40%

comparaison03 <- data.frame(matrix(ncol = 5, 
                                   nrow = 42))

colnames(comparaison03) <- c(
  "Estimation", 
  "Performance",
  "Scénario06_50%",
  "Scénario06_70%",
  "Scénario06_85%"
)

comparaison03$Estimation <- c("RegLog", "-", "-", "-", "-","-", "IPW", "-", "-", "-", "-", "-", 
                               "TNDDR_RF", "-", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-", "-",
                               "TNDDR_Mars", "-", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-", "-",
                               "TNDDR_GLM", "-", "-", "-", "-", "-")

comparaison03$Performance <- rep(c("Biais_med", "Biais_moy", "Variance", "MSE", "Précision", "%Cov"), 7)

comparaison03$`Scénario06_50%` <- c(
  
  (median(resultats_6_1$RRc) - mean(l_vraiRRc_6_1)), Tab01_6_1$Autres, Tab01_6_1$`Erreur de Monte Carlo`[4],
  (median(resultats2_6_1$RRm) - mean(l_vraiRRm_6_1)), Tab02_6_1$Autres, Tab02_6_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_1$RRm_RF) - mean(l_vraiRRm_6_1)), Tab03_6_1$Autres[1:4], Tab03_6_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_1$RRm_Lasso) - mean(l_vraiRRm_6_1)),Tab03_6_1$Autres[5:8], Tab03_6_1$`Erreur de Monte Carlo`[8],
  (median(resultats3_6_1$RRm_Mars) - mean(l_vraiRRm_6_1)), Tab03_6_1$Autres[9:12], Tab03_6_1$`Erreur de Monte Carlo`[12],
  (median(resultats3_6_1$RRm_RN) - mean(l_vraiRRm_6_1)), Tab03_6_1$Autres[13:16], Tab03_6_1$`Erreur de Monte Carlo`[16],
  (median(resultats3_6_1$RRm_GLM) - mean(l_vraiRRm_6_1)), Tab03_6_1$Autres[17:20], Tab03_6_1$`Erreur de Monte Carlo`[20]
  
)

comparaison03$`Scénario06_70%` <- c(
  
  (median(resultats_6_2$RRc) - mean(l_vraiRRc_6_2)), Tab01_6_2$Autres, Tab01_6_2$`Erreur de Monte Carlo`[4],
  (median(resultats2_6_2$RRm) - mean(l_vraiRRm_6_2)), Tab02_6_2$Autres, Tab02_6_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_2$RRm_RF) - mean(l_vraiRRm_6_2)), Tab03_6_2$Autres[1:4], Tab03_6_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_2$RRm_Lasso) - mean(l_vraiRRm_6_2)),Tab03_6_2$Autres[5:8], Tab03_6_2$`Erreur de Monte Carlo`[8],
  (median(resultats3_6_2$RRm_Mars) - mean(l_vraiRRm_6_2)), Tab03_6_2$Autres[9:12], Tab03_6_2$`Erreur de Monte Carlo`[12],
  (median(resultats3_6_2$RRm_RN) - mean(l_vraiRRm_6_2)), Tab03_6_2$Autres[13:16], Tab03_6_2$`Erreur de Monte Carlo`[16],
  (median(resultats3_6_2$RRm_GLM) - mean(l_vraiRRm_6_2)), Tab03_6_2$Autres[17:20], Tab03_6_2$`Erreur de Monte Carlo`[20]
  
)

comparaison03$`Scénario06_85%` <- c(
  
  (median(resultats_6_3$RRc) - mean(l_vraiRRc_6_3)), Tab01_6_3$Autres, Tab01_6_3$`Erreur de Monte Carlo`[4],
  (median(resultats2_6_3$RRm) - mean(l_vraiRRm_6_3)), Tab02_6_3$Autres, Tab02_6_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_3$RRm_RF) - mean(l_vraiRRm_6_3)), Tab03_6_3$Autres[1:4], Tab03_6_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_3$RRm_Lasso, na.rm = TRUE) - mean(l_vraiRRm_6_3)),Tab03_6_3$Autres[5:8], Tab03_6_3$`Erreur de Monte Carlo`[8],
  (median(resultats3_6_3$RRm_Mars) - mean(l_vraiRRm_6_3)), Tab03_6_3$Autres[9:12], Tab03_6_3$`Erreur de Monte Carlo`[12],
  (median(resultats3_6_3$RRm_RN) - mean(l_vraiRRm_6_3)), Tab03_6_3$Autres[13:16], Tab03_6_3$`Erreur de Monte Carlo`[16],
  (median(resultats3_6_3$RRm_GLM) - mean(l_vraiRRm_6_3)), Tab03_6_3$Autres[17:20], Tab03_6_3$`Erreur de Monte Carlo`[20]
  
)

comparaison03$`Scénario06_50%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison03$`Scénario06_50%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison03$`Scénario06_70%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison03$`Scénario06_70%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison03$`Scénario06_85%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison03$`Scénario06_85%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison03)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison03$`Scénario06_50%`[[2]], comparaison03$`Scénario06_70%`[[2]], comparaison03$`Scénario06_85%`[[2]],
                               comparaison03$`Scénario06_50%`[[8]], comparaison03$`Scénario06_70%`[[8]], comparaison03$`Scénario06_85%`[[8]],
                               comparaison03$`Scénario06_50%`[[14]], comparaison03$`Scénario06_70%`[[14]], comparaison03$`Scénario06_85%`[[14]],
                               comparaison03$`Scénario06_50%`[[20]], comparaison03$`Scénario06_70%`[[20]], comparaison03$`Scénario06_85%`[[20]],
                               comparaison03$`Scénario06_50%`[[26]], comparaison03$`Scénario06_70%`[[26]], comparaison03$`Scénario06_85%`[[26]],
                               comparaison03$`Scénario06_50%`[[32]], comparaison03$`Scénario06_70%`[[32]], comparaison03$`Scénario06_85%`[[32]],
                               comparaison03$`Scénario06_50%`[[38]], comparaison03$`Scénario06_70%`[[38]], comparaison03$`Scénario06_85%`[[38]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                   Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture

dat <- data.frame(x_cov = c(comparaison03$`Scénario06_50%`[[6]], comparaison03$`Scénario06_70%`[[6]], comparaison03$`Scénario06_85%`[[6]],
                               comparaison03$`Scénario06_50%`[[12]], comparaison03$`Scénario06_70%`[[12]], comparaison03$`Scénario06_85%`[[12]],
                               comparaison03$`Scénario06_50%`[[18]], comparaison03$`Scénario06_70%`[[18]], comparaison03$`Scénario06_85%`[[18]],
                               comparaison03$`Scénario06_50%`[[24]], comparaison03$`Scénario06_70%`[[24]], comparaison03$`Scénario06_85%`[[24]],
                               comparaison03$`Scénario06_50%`[[30]], comparaison03$`Scénario06_70%`[[30]], comparaison03$`Scénario06_85%`[[30]],
                               comparaison03$`Scénario06_50%`[[36]], comparaison03$`Scénario06_70%`[[36]], comparaison03$`Scénario06_85%`[[36]],
                               comparaison03$`Scénario06_50%`[[42]], comparaison03$`Scénario06_70%`[[42]], comparaison03$`Scénario06_85%`[[42]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                   Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Biais median

dat <- data.frame(x_biais = c(comparaison03$`Scénario06_50%`[[1]], comparaison03$`Scénario06_70%`[[1]], comparaison03$`Scénario06_85%`[[1]],
                              comparaison03$`Scénario06_50%`[[7]], comparaison03$`Scénario06_70%`[[7]], comparaison03$`Scénario06_85%`[[7]],
                              comparaison03$`Scénario06_50%`[[13]], comparaison03$`Scénario06_70%`[[13]], comparaison03$`Scénario06_85%`[[13]],
                              comparaison03$`Scénario06_50%`[[19]], comparaison03$`Scénario06_70%`[[19]], comparaison03$`Scénario06_85%`[[19]],
                              comparaison03$`Scénario06_50%`[[25]], comparaison03$`Scénario06_70%`[[25]], comparaison03$`Scénario06_85%`[[25]],
                              comparaison03$`Scénario06_50%`[[31]], comparaison03$`Scénario06_70%`[[31]], comparaison03$`Scénario06_85%`[[31]],
                              comparaison03$`Scénario06_50%`[[37]], comparaison03$`Scénario06_70%`[[37]], comparaison03$`Scénario06_85%`[[37]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

################################################################################
################################################################################
# Scénarios avec différentes valeurs de prévalence de I1 et un taux de co-infection entre 30% et 40%

## Chargement des résultats bruts

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario07_10%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario07_30%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario07_50%.RData")

# nsim = 1000, Scénarios 07 prévalence de I1  ~ 10%, 30% et 50%, taux de co-infection ~30%-40%

comparaison04 <- data.frame(matrix(ncol = 5, 
                                   nrow = 42))

colnames(comparaison04) <- c(
  "Estimation", 
  "Performance",
  "Scénario07_10%",
  "Scénario07_20%",
  "Scénario07_50%"
)

comparaison04$Estimation <- c("RegLog", "-", "-", "-", "-","-", "IPW", "-", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-", "-",
                              "TNDDR_GLM", "-", "-", "-", "-", "-")

comparaison04$Performance <- rep(c("Biais_med", "Biais_moy", "Variance", "MSE", "Précision", "%Cov"), 7)

comparaison04$`Scénario07_10%` <- c(
  
  (median(resultats_7_1$RRc) - mean(l_vraiRRc_7_1)), Tab01_7_1$Autres, Tab01_7_1$`Erreur de Monte Carlo`[4],
  (median(resultats2_7_1$RRm) - mean(l_vraiRRm_7_1)), Tab02_7_1$Autres, Tab02_7_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_1$RRm_RF) - mean(l_vraiRRm_7_1)), Tab03_7_1$Autres[1:4], Tab03_7_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_1$RRm_Lasso, na.rm = TRUE) - mean(l_vraiRRm_7_1)),Tab03_7_1$Autres[5:8], Tab03_7_1$`Erreur de Monte Carlo`[8],
  (median(resultats3_7_1$RRm_Mars) - mean(l_vraiRRm_7_1)), Tab03_7_1$Autres[9:12], Tab03_7_1$`Erreur de Monte Carlo`[12],
  (median(resultats3_7_1$RRm_RN) - mean(l_vraiRRm_7_1)), Tab03_7_1$Autres[13:16], Tab03_7_1$`Erreur de Monte Carlo`[16],
  (median(resultats3_7_1$RRm_GLM) - mean(l_vraiRRm_7_1)), Tab03_7_1$Autres[17:20], Tab03_7_1$`Erreur de Monte Carlo`[20]
  
)

comparaison04$`Scénario07_20%` <- c(
  
  (median(resultats_7_2$RRc) - mean(l_vraiRRc_7_2)), Tab01_7_2$Autres, Tab01_7_2$`Erreur de Monte Carlo`[4],
  (median(resultats2_7_2$RRm) - mean(l_vraiRRm_7_2)), Tab02_7_2$Autres, Tab02_7_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_2$RRm_RF) - mean(l_vraiRRm_7_2)), Tab03_7_2$Autres[1:4], Tab03_7_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_2$RRm_Lasso) - mean(l_vraiRRm_7_2)),Tab03_7_2$Autres[5:8], Tab03_7_2$`Erreur de Monte Carlo`[8],
  (median(resultats3_7_2$RRm_Mars) - mean(l_vraiRRm_7_2)), Tab03_7_2$Autres[9:12], Tab03_7_2$`Erreur de Monte Carlo`[12],
  (median(resultats3_7_2$RRm_RN) - mean(l_vraiRRm_7_2)), Tab03_7_2$Autres[13:16], Tab03_7_2$`Erreur de Monte Carlo`[16],
  (median(resultats3_7_2$RRm_GLM) - mean(l_vraiRRm_7_2)), Tab03_7_2$Autres[17:20], Tab03_7_2$`Erreur de Monte Carlo`[20]
  
)

comparaison04$`Scénario07_50%` <- c(
  
  (median(resultats_7_3$RRc) - mean(l_vraiRRc_7_3)), Tab01_7_3$Autres, Tab01_7_3$`Erreur de Monte Carlo`[4],
  (median(resultats2_7_3$RRm) - mean(l_vraiRRm_7_3)), Tab02_7_3$Autres, Tab02_7_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_3$RRm_RF) - mean(l_vraiRRm_7_3)), Tab03_7_3$Autres[1:4], Tab03_7_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_3$RRm_Lasso) - mean(l_vraiRRm_7_3)),Tab03_7_3$Autres[5:8], Tab03_7_3$`Erreur de Monte Carlo`[8],
  (median(resultats3_7_3$RRm_Mars) - mean(l_vraiRRm_7_3)), Tab03_7_3$Autres[9:12], Tab03_7_3$`Erreur de Monte Carlo`[12],
  (median(resultats3_7_3$RRm_RN) - mean(l_vraiRRm_7_3)), Tab03_7_3$Autres[13:16], Tab03_7_3$`Erreur de Monte Carlo`[16],
  (median(resultats3_7_3$RRm_GLM) - mean(l_vraiRRm_7_3)), Tab03_7_3$Autres[17:20], Tab03_7_3$`Erreur de Monte Carlo`[20]
  
)

comparaison04$`Scénario07_10%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison04$`Scénario07_10%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison04$`Scénario07_20%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison04$`Scénario07_20%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison04$`Scénario07_50%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison04$`Scénario07_50%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison04)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison04$`Scénario07_10%`[[2]], comparaison04$`Scénario07_20%`[[2]], comparaison04$`Scénario07_50%`[[2]],
                               comparaison04$`Scénario07_10%`[[8]], comparaison04$`Scénario07_20%`[[8]], comparaison04$`Scénario07_50%`[[8]],
                               comparaison04$`Scénario07_10%`[[14]], comparaison04$`Scénario07_20%`[[14]], comparaison04$`Scénario07_50%`[[14]],
                               comparaison04$`Scénario07_10%`[[20]], comparaison04$`Scénario07_20%`[[20]], comparaison04$`Scénario07_50%`[[20]],
                               comparaison04$`Scénario07_10%`[[26]], comparaison04$`Scénario07_20%`[[26]], comparaison04$`Scénario07_50%`[[26]],
                               comparaison04$`Scénario07_10%`[[32]], comparaison04$`Scénario07_20%`[[32]], comparaison04$`Scénario07_50%`[[32]],
                               comparaison04$`Scénario07_10%`[[38]], comparaison04$`Scénario07_20%`[[38]], comparaison04$`Scénario07_50%`[[38]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                   Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture

dat <- data.frame(x_cov = c(comparaison04$`Scénario07_10%`[[6]], comparaison04$`Scénario07_20%`[[6]], comparaison04$`Scénario07_50%`[[6]],
                               comparaison04$`Scénario07_10%`[[12]], comparaison04$`Scénario07_20%`[[12]], comparaison04$`Scénario07_50%`[[12]],
                               comparaison04$`Scénario07_10%`[[18]], comparaison04$`Scénario07_20%`[[18]], comparaison04$`Scénario07_50%`[[18]],
                               comparaison04$`Scénario07_10%`[[24]], comparaison04$`Scénario07_20%`[[24]], comparaison04$`Scénario07_50%`[[24]],
                               comparaison04$`Scénario07_10%`[[30]], comparaison04$`Scénario07_20%`[[30]], comparaison04$`Scénario07_50%`[[30]],
                               comparaison04$`Scénario07_10%`[[36]], comparaison04$`Scénario07_20%`[[36]], comparaison04$`Scénario07_50%`[[36]],
                               comparaison04$`Scénario07_10%`[[42]], comparaison04$`Scénario07_20%`[[42]], comparaison04$`Scénario07_50%`[[42]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                   Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Biais median

dat <- data.frame(x_biais = c(comparaison04$`Scénario07_10%`[[1]], comparaison04$`Scénario07_20%`[[1]], comparaison04$`Scénario07_50%`[[1]],
                              comparaison04$`Scénario07_10%`[[7]], comparaison04$`Scénario07_20%`[[7]], comparaison04$`Scénario07_50%`[[7]],
                              comparaison04$`Scénario07_10%`[[13]], comparaison04$`Scénario07_20%`[[13]], comparaison04$`Scénario07_50%`[[13]],
                              comparaison04$`Scénario07_10%`[[19]], comparaison04$`Scénario07_20%`[[19]], comparaison04$`Scénario07_50%`[[19]],
                              comparaison04$`Scénario07_10%`[[25]], comparaison04$`Scénario07_20%`[[25]], comparaison04$`Scénario07_50%`[[25]],
                              comparaison04$`Scénario07_10%`[[31]], comparaison04$`Scénario07_20%`[[31]], comparaison04$`Scénario07_50%`[[31]],
                              comparaison04$`Scénario07_10%`[[37]], comparaison04$`Scénario07_20%`[[37]], comparaison04$`Scénario07_50%`[[37]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

################################################################################
################################################################################
# Scénarios avec différentes valeurs de prévalence de W2 et un taux de co-infection entre 30% et 40%

## Chargement des résultats bruts

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario08_5%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario08_10%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario08_15%.RData")

# nsim = 1000, Scénarios 08 prévalence de W2  ~ 5%, 10% et 15%, taux de co-infection ~30%-40%

# Tableau de comparaison

comparaison05 <- data.frame(matrix(ncol = 5, 
                                   nrow = 42))

colnames(comparaison05) <- c(
  "Estimation", 
  "Performance",
  "Scénario08_5%",
  "Scénario08_10%",
  "Scénario08_15%"
)

comparaison05$Estimation <- c("RegLog", "-", "-", "-", "-","-", "IPW", "-", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-", "-",
                              "TNDDR_GLM", "-", "-", "-", "-", "-")

comparaison05$Performance <- rep(c("Biais_med", "Biais_moy", "Variance", "MSE", "Précision", "%Cov"), 7)

comparaison05$`Scénario08_5%` <- c(
  
  (median(resultats_8_1$RRc) - mean(l_vraiRRc_8_1)), Tab01_8_1$Autres, Tab01_8_1$`Erreur de Monte Carlo`[4],
  (median(resultats2_8_1$RRm) - mean(l_vraiRRm_8_1)), Tab02_8_1$Autres, Tab02_8_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_1$RRm_RF) - mean(l_vraiRRm_8_1)), Tab03_8_1$Autres[1:4], Tab03_8_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_1$RRm_Lasso) - mean(l_vraiRRm_8_1)),Tab03_8_1$Autres[5:8], Tab03_8_1$`Erreur de Monte Carlo`[8],
  (median(resultats3_8_1$RRm_Mars) - mean(l_vraiRRm_8_1)), Tab03_8_1$Autres[9:12], Tab03_8_1$`Erreur de Monte Carlo`[12],
  (median(resultats3_8_1$RRm_RN) - mean(l_vraiRRm_8_1)), Tab03_8_1$Autres[13:16], Tab03_8_1$`Erreur de Monte Carlo`[16],
  (median(resultats3_8_1$RRm_GLM) - mean(l_vraiRRm_8_1)), Tab03_8_1$Autres[17:20], Tab03_8_1$`Erreur de Monte Carlo`[20]
  
)

comparaison05$`Scénario08_10%` <- c(
  
  (median(resultats_8_2$RRc) - mean(l_vraiRRc_8_2)), Tab01_8_2$Autres, Tab01_8_2$`Erreur de Monte Carlo`[4],
  (median(resultats2_8_2$RRm) - mean(l_vraiRRm_8_2)), Tab02_8_2$Autres, Tab02_8_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_2$RRm_RF) - mean(l_vraiRRm_8_2)), Tab03_8_2$Autres[1:4], Tab03_8_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_2$RRm_Lasso) - mean(l_vraiRRm_8_2)),Tab03_8_2$Autres[5:8], Tab03_8_2$`Erreur de Monte Carlo`[8],
  (median(resultats3_8_2$RRm_Mars) - mean(l_vraiRRm_8_2)), Tab03_8_2$Autres[9:12], Tab03_8_2$`Erreur de Monte Carlo`[12],
  (median(resultats3_8_2$RRm_RN) - mean(l_vraiRRm_8_2)), Tab03_8_2$Autres[13:16], Tab03_8_2$`Erreur de Monte Carlo`[16],
  (median(resultats3_8_2$RRm_GLM) - mean(l_vraiRRm_8_2)), Tab03_8_2$Autres[17:20], Tab03_8_2$`Erreur de Monte Carlo`[20]
  
)

comparaison05$`Scénario08_15%` <- c(
  
  (median(resultats_8_3$RRc) - mean(l_vraiRRc_8_3)), Tab01_8_3$Autres, Tab01_8_3$`Erreur de Monte Carlo`[4],
  (median(resultats2_8_3$RRm) - mean(l_vraiRRm_8_3)), Tab02_8_3$Autres, Tab02_8_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_3$RRm_RF) - mean(l_vraiRRm_8_3)), Tab03_8_3$Autres[1:4], Tab03_8_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_3$RRm_Lasso) - mean(l_vraiRRm_8_3)),Tab03_8_3$Autres[5:8], Tab03_8_3$`Erreur de Monte Carlo`[8],
  (median(resultats3_8_3$RRm_Mars) - mean(l_vraiRRm_8_3)), Tab03_8_3$Autres[9:12], Tab03_8_3$`Erreur de Monte Carlo`[12],
  (median(resultats3_8_3$RRm_RN) - mean(l_vraiRRm_8_3)), Tab03_8_3$Autres[13:16], Tab03_8_3$`Erreur de Monte Carlo`[16],
  (median(resultats3_8_3$RRm_GLM) - mean(l_vraiRRm_8_3)), Tab03_8_3$Autres[17:20], Tab03_8_3$`Erreur de Monte Carlo`[20]
  
)

comparaison05$`Scénario08_5%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison05$`Scénario08_5%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison05$`Scénario08_10%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison05$`Scénario08_10%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison05$`Scénario08_15%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison05$`Scénario08_15%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison05)

####################### Représentation graphique ###############################

# Biais moyen 

dat <- data.frame(x_biais = c(comparaison05$`Scénario08_5%`[[2]], comparaison05$`Scénario08_10%`[[2]], comparaison05$`Scénario08_15%`[[2]],
                              comparaison05$`Scénario08_5%`[[8]], comparaison05$`Scénario08_10%`[[8]], comparaison05$`Scénario08_15%`[[8]],
                              comparaison05$`Scénario08_5%`[[14]], comparaison05$`Scénario08_10%`[[14]], comparaison05$`Scénario08_15%`[[14]],
                              comparaison05$`Scénario08_5%`[[20]], comparaison05$`Scénario08_10%`[[20]], comparaison05$`Scénario08_15%`[[20]],
                              comparaison05$`Scénario08_5%`[[26]], comparaison05$`Scénario08_10%`[[26]], comparaison05$`Scénario08_15%`[[26]],
                              comparaison05$`Scénario08_5%`[[32]], comparaison05$`Scénario08_10%`[[32]], comparaison05$`Scénario08_15%`[[32]],
                              comparaison05$`Scénario08_5%`[[38]], comparaison05$`Scénario08_10%`[[38]], comparaison05$`Scénario08_15%`[[38]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture

dat <- data.frame(x_cov = c(comparaison05$`Scénario08_10%`[[6]], comparaison05$`Scénario08_15%`[[6]], comparaison05$`Scénario08_20%`[[6]],
                            comparaison05$`Scénario08_10%`[[12]], comparaison05$`Scénario08_15%`[[12]], comparaison05$`Scénario08_20%`[[12]],
                            comparaison05$`Scénario08_10%`[[18]], comparaison05$`Scénario08_15%`[[18]], comparaison05$`Scénario08_20%`[[18]],
                            comparaison05$`Scénario08_10%`[[24]], comparaison05$`Scénario08_15%`[[24]], comparaison05$`Scénario08_20%`[[24]],
                            comparaison05$`Scénario08_10%`[[30]], comparaison05$`Scénario08_15%`[[30]], comparaison05$`Scénario08_20%`[[30]],
                            comparaison05$`Scénario08_10%`[[36]], comparaison05$`Scénario08_15%`[[36]], comparaison05$`Scénario08_20%`[[36]],
                            comparaison05$`Scénario08_10%`[[42]], comparaison05$`Scénario08_15%`[[42]], comparaison05$`Scénario08_20%`[[42]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Biais median

dat <- data.frame(x_biais = c(comparaison05$`Scénario08_5%`[[1]], comparaison05$`Scénario08_10%`[[1]], comparaison05$`Scénario08_15%`[[1]],
                              comparaison05$`Scénario08_5%`[[7]], comparaison05$`Scénario08_10%`[[7]], comparaison05$`Scénario08_15%`[[7]],
                              comparaison05$`Scénario08_5%`[[13]], comparaison05$`Scénario08_10%`[[13]], comparaison05$`Scénario08_15%`[[13]],
                              comparaison05$`Scénario08_5%`[[19]], comparaison05$`Scénario08_10%`[[19]], comparaison05$`Scénario08_15%`[[19]],
                              comparaison05$`Scénario08_5%`[[25]], comparaison05$`Scénario08_10%`[[25]], comparaison05$`Scénario08_15%`[[25]],
                              comparaison05$`Scénario08_5%`[[31]], comparaison05$`Scénario08_10%`[[31]], comparaison05$`Scénario08_15%`[[31]],
                              comparaison05$`Scénario08_5%`[[37]], comparaison05$`Scénario08_10%`[[37]], comparaison05$`Scénario08_15%`[[37]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

################################################################################
################################################################################
# Scénarios avec différentes valeurs de prévalence de W1 et un taux de co-infection entre 30% et 40%

## Chargement des résultats bruts

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario09_10%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario09_15%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario09_20%.RData")

# nsim = 1000, Scénarios 09 prévalence de W1  ~ 10%, 15% et 30%, taux de co-infection ~30%-40%

# Tableau de comparaison

comparaison06 <- data.frame(matrix(ncol = 5, 
                                   nrow = 42))

colnames(comparaison06) <- c(
  "Estimation", 
  "Performance",
  "Scénario09_10%",
  "Scénario09_15%",
  "Scénario09_20%"
)

comparaison06$Estimation <- c("RegLog", "-", "-", "-", "-","-", "IPW", "-", "-", "-", "-", "-", 
                               "TNDDR_RF", "-", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-", "-",
                               "TNDDR_Mars", "-", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-", "-",
                               "TNDDR_GLM", "-", "-", "-", "-", "-")

comparaison06$Performance <- rep(c("Biais_med", "Biais_moy", "Variance", "MSE", "Précision", "%Cov"), 7)

comparaison06$`Scénario09_10%` <- c(
  
  (median(resultats_9_1$RRc) - mean(l_vraiRRc_9_1)), Tab01_9_1$Autres, Tab01_9_1$`Erreur de Monte Carlo`[4],
  (median(resultats2_9_1$RRm) - mean(l_vraiRRm_9_1)), Tab02_9_1$Autres, Tab02_9_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_1$RRm_RF) - mean(l_vraiRRm_9_1)), Tab03_9_1$Autres[1:4], Tab03_9_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_1$RRm_Lasso) - mean(l_vraiRRm_9_1)),Tab03_9_1$Autres[5:8], Tab03_9_1$`Erreur de Monte Carlo`[8],
  (median(resultats3_9_1$RRm_Mars) - mean(l_vraiRRm_9_1)), Tab03_9_1$Autres[9:12], Tab03_9_1$`Erreur de Monte Carlo`[12],
  (median(resultats3_9_1$RRm_RN) - mean(l_vraiRRm_9_1)), Tab03_9_1$Autres[13:16], Tab03_9_1$`Erreur de Monte Carlo`[16],
  (median(resultats3_9_1$RRm_GLM) - mean(l_vraiRRm_9_1)), Tab03_9_1$Autres[17:20], Tab03_9_1$`Erreur de Monte Carlo`[20]
  
)

comparaison06$`Scénario09_15%` <- c(
  
  (median(resultats_9_2$RRc) - mean(l_vraiRRc_9_2)), Tab01_9_2$Autres, Tab01_9_2$`Erreur de Monte Carlo`[4],
  (median(resultats2_9_2$RRm) - mean(l_vraiRRm_9_2)), Tab02_9_2$Autres, Tab02_9_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_2$RRm_RF) - mean(l_vraiRRm_9_2)), Tab03_9_2$Autres[1:4], Tab03_9_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_2$RRm_Lasso) - mean(l_vraiRRm_9_2)),Tab03_9_2$Autres[5:8], Tab03_9_2$`Erreur de Monte Carlo`[8],
  (median(resultats3_9_2$RRm_Mars) - mean(l_vraiRRm_9_2)), Tab03_9_2$Autres[9:12], Tab03_9_2$`Erreur de Monte Carlo`[12],
  (median(resultats3_9_2$RRm_RN) - mean(l_vraiRRm_9_2)), Tab03_9_2$Autres[13:16], Tab03_9_2$`Erreur de Monte Carlo`[16],
  (median(resultats3_9_2$RRm_GLM) - mean(l_vraiRRm_9_2)), Tab03_9_2$Autres[17:20], Tab03_9_2$`Erreur de Monte Carlo`[20]
  
)

comparaison06$`Scénario09_20%` <- c(
  
  (median(resultats_9_3$RRc) - mean(l_vraiRRc_9_3)), Tab01_9_3$Autres, Tab01_9_3$`Erreur de Monte Carlo`[4],
  (median(resultats2_9_3$RRm) - mean(l_vraiRRm_9_3)), Tab02_9_3$Autres, Tab02_9_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_3$RRm_RF) - mean(l_vraiRRm_9_3)), Tab03_9_3$Autres[1:4], Tab03_9_3$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_3$RRm_Lasso) - mean(l_vraiRRm_9_3)),Tab03_9_3$Autres[5:8], Tab03_9_3$`Erreur de Monte Carlo`[8],
  (median(resultats3_9_3$RRm_Mars) - mean(l_vraiRRm_9_3)), Tab03_9_3$Autres[9:12], Tab03_9_3$`Erreur de Monte Carlo`[12],
  (median(resultats3_9_3$RRm_RN) - mean(l_vraiRRm_9_3)), Tab03_9_3$Autres[13:16], Tab03_9_3$`Erreur de Monte Carlo`[16],
  (median(resultats3_9_3$RRm_GLM) - mean(l_vraiRRm_9_3)), Tab03_9_3$Autres[17:20], Tab03_9_3$`Erreur de Monte Carlo`[20]
  
)

comparaison06$`Scénario09_10%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison06$`Scénario09_10%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison06$`Scénario09_15%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison06$`Scénario09_15%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison06$`Scénario09_20%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison06$`Scénario09_20%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison06)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison06$`Scénario09_10%`[[2]], comparaison06$`Scénario09_15%`[[2]], comparaison06$`Scénario09_20%`[[2]],
                              comparaison06$`Scénario09_10%`[[8]], comparaison06$`Scénario09_15%`[[8]], comparaison06$`Scénario09_20%`[[8]],
                              comparaison06$`Scénario09_10%`[[14]], comparaison06$`Scénario09_15%`[[14]], comparaison06$`Scénario09_20%`[[14]],
                              comparaison06$`Scénario09_10%`[[20]], comparaison06$`Scénario09_15%`[[20]], comparaison06$`Scénario09_20%`[[20]],
                              comparaison06$`Scénario09_10%`[[26]], comparaison06$`Scénario09_15%`[[26]], comparaison06$`Scénario09_20%`[[26]],
                              comparaison06$`Scénario09_10%`[[32]], comparaison06$`Scénario09_15%`[[30]], comparaison06$`Scénario09_20%`[[30]],
                              comparaison06$`Scénario09_10%`[[38]], comparaison06$`Scénario09_15%`[[38]], comparaison06$`Scénario09_20%`[[38]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture

dat <- data.frame(x_cov = c(comparaison06$`Scénario09_10%`[[6]], comparaison06$`Scénario09_15%`[[6]], comparaison06$`Scénario09_20%`[[6]],
                            comparaison06$`Scénario09_10%`[[12]], comparaison06$`Scénario09_15%`[[12]], comparaison06$`Scénario09_20%`[[12]],
                            comparaison06$`Scénario09_10%`[[18]], comparaison06$`Scénario09_15%`[[18]], comparaison06$`Scénario09_20%`[[18]],
                            comparaison06$`Scénario09_10%`[[24]], comparaison06$`Scénario09_15%`[[24]], comparaison06$`Scénario09_20%`[[24]],
                            comparaison06$`Scénario09_10%`[[30]], comparaison06$`Scénario09_15%`[[30]], comparaison06$`Scénario09_20%`[[30]],
                            comparaison06$`Scénario09_10%`[[36]], comparaison06$`Scénario09_15%`[[36]], comparaison06$`Scénario09_20%`[[36]],
                            comparaison06$`Scénario09_10%`[[42]], comparaison06$`Scénario09_15%`[[42]], comparaison06$`Scénario09_20%`[[42]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Biais median

dat <- data.frame(x_biais = c(comparaison06$`Scénario09_10%`[[1]], comparaison06$`Scénario09_15%`[[1]], comparaison06$`Scénario09_20%`[[1]],
                              comparaison06$`Scénario09_10%`[[7]], comparaison06$`Scénario09_15%`[[7]], comparaison06$`Scénario09_20%`[[7]],
                              comparaison06$`Scénario09_10%`[[13]], comparaison06$`Scénario09_15%`[[13]], comparaison06$`Scénario09_20%`[[13]],
                              comparaison06$`Scénario09_10%`[[19]], comparaison06$`Scénario09_15%`[[19]], comparaison06$`Scénario09_20%`[[19]],
                              comparaison06$`Scénario09_10%`[[25]], comparaison06$`Scénario09_15%`[[25]], comparaison06$`Scénario09_20%`[[25]],
                              comparaison06$`Scénario09_10%`[[31]], comparaison06$`Scénario09_15%`[[31]], comparaison06$`Scénario09_20%`[[31]],
                              comparaison06$`Scénario09_10%`[[37]], comparaison06$`Scénario09_15%`[[37]], comparaison06$`Scénario09_20%`[[37]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

################################################################################
################################################################################
## Taille d'échantillo = 5000
################################################################################
## Scénarios avec différents taux de co-infection : 0%, 10%, 20% et 40%
# nsim = 1000, Scénarios 01, 02, 03 et 04

# Chargement des résultats bruts

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario01_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario02_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario03_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats bruts/Résultats_Scénario04_5000.RData")

# Tableau de comparaison 

comparaison01 <- data.frame(matrix(ncol = 6, 
                                   nrow = 42))

colnames(comparaison01) <- c(
  "Estimation", 
  "Performance",
  "Scénario01_0%",
  "Scénario02_10%",
  "Scénario03_20%",
  "Scénario04_40%"
)

comparaison01$Estimation <- c("RegLog", "-", "-", "-", "-","-", "IPW", "-", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-", "-",
                              "TNDDR_GLM", "-", "-", "-", "-", "-")

comparaison01$Performance <- rep(c("Biais_med", "Biais_moy", "Variance", "MSE", "Précision", "%Cov"), 7)

comparaison01$`Scénario01_0%` <- c(
  
  (median(resultats_1_5000$RRc) - mean(l_vraiRRc_1_5000)), Tab01_1_5000$Autres, Tab01_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_1_5000$RRm) - mean(l_vraiRRm_1_5000)), Tab02_1_5000$Autres, Tab02_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_1_5000$RRm_RF) - mean(l_vraiRRm_1_5000)), Tab03_1_5000$Autres[1:4], Tab03_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_1_5000$RRm_Lasso) - mean(l_vraiRRm_1_5000)),Tab03_1_5000$Autres[5:8], Tab03_1_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_1_5000$RRm_Mars) - mean(l_vraiRRm_1_5000)), Tab03_1_5000$Autres[9:12], Tab03_1_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_1_5000$RRm_RN) - mean(l_vraiRRm_1_5000)), Tab03_1_5000$Autres[13:16], Tab03_1_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_1_5000$RRm_GLM) - mean(l_vraiRRm_1_5000)), Tab03_1_5000$Autres[17:20], Tab03_1_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario02_10%` <- c(
  
  (median(resultats_2_5000$RRc) - mean(l_vraiRRc_2_5000)), Tab01_2_5000$Autres, Tab01_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_2_5000$RRm) - mean(l_vraiRRm_2_5000)), Tab02_2_5000$Autres, Tab02_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_2_5000$RRm_RF) - mean(l_vraiRRm_2_5000)), Tab03_2_5000$Autres[1:4], Tab03_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_2_5000$RRm_Lasso) - mean(l_vraiRRm_2_5000)),Tab03_2_5000$Autres[5:8], Tab03_2_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_2_5000$RRm_Mars) - mean(l_vraiRRm_2_5000)), Tab03_2_5000$Autres[9:12], Tab03_2_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_2_5000$RRm_RN) - mean(l_vraiRRm_2_5000)), Tab03_2_5000$Autres[13:16], Tab03_2_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_2_5000$RRm_GLM) - mean(l_vraiRRm_2_5000)), Tab03_2_5000$Autres[17:20], Tab03_2_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario03_20%` <- c(
  
  (median(resultats_3_5000$RRc) - mean(l_vraiRRc_3_5000)), Tab01_3_5000$Autres, Tab01_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_3_5000$RRm) - mean(l_vraiRRm_3_5000)), Tab02_3_5000$Autres, Tab02_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_3_5000$RRm_RF) - mean(l_vraiRRm_3_5000)), Tab03_3_5000$Autres[1:4], Tab03_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_3_5000$RRm_Lasso) - mean(l_vraiRRm_3_5000)),Tab03_3_5000$Autres[5:8], Tab03_3_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_3_5000$RRm_Mars) - mean(l_vraiRRm_3_5000)), Tab03_3_5000$Autres[9:12], Tab03_3_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_3_5000$RRm_RN) - mean(l_vraiRRm_3_5000)), Tab03_3_5000$Autres[13:16], Tab03_3_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_3_5000$RRm_GLM) - mean(l_vraiRRm_3_5000)), Tab03_3_5000$Autres[17:20], Tab03_3_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario04_40%` <- c(
  
  (median(resultats_4_5000$RRc) - mean(l_vraiRRc_4_5000)), Tab01_4_5000$Autres, Tab01_4_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_4_5000$RRm) - mean(l_vraiRRm_4_5000)), Tab02_4_5000$Autres, Tab02_4_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_4_5000$RRm_RF) - mean(l_vraiRRm_4_5000)), Tab03_4_5000$Autres[1:4], Tab03_4_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_4_5000$RRm_Lasso) - mean(l_vraiRRm_4_5000)),Tab03_4_5000$Autres[5:8], Tab03_4_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_4_5000$RRm_Mars) - mean(l_vraiRRm_4_5000)), Tab03_4_5000$Autres[9:12], Tab03_4_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_4_5000$RRm_RN) - mean(l_vraiRRm_4_5000)), Tab03_4_5000$Autres[13:16], Tab03_4_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_4_5000$RRm_GLM) - mean(l_vraiRRm_4_5000)), Tab03_4_5000$Autres[17:20], Tab03_4_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario01_0%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario01_0%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison01$`Scénario02_10%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario02_10%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison01$`Scénario03_20%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario03_20%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison01$`Scénario04_40%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario04_40%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison01)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison01$`Scénario01_0%`[[2]], comparaison01$`Scénario02_10%`[[2]], comparaison01$`Scénario03_20%`[[2]], comparaison01$`Scénario04_40%`[[2]],
                              comparaison01$`Scénario01_0%`[[8]], comparaison01$`Scénario02_10%`[[8]], comparaison01$`Scénario03_20%`[[8]], comparaison01$`Scénario04_40%`[[8]],
                              comparaison01$`Scénario01_0%`[[14]], comparaison01$`Scénario02_10%`[[14]], comparaison01$`Scénario03_20%`[[14]], comparaison01$`Scénario04_40%`[[14]],
                              comparaison01$`Scénario01_0%`[[20]], comparaison01$`Scénario02_10%`[[20]], comparaison01$`Scénario03_20%`[[20]], comparaison01$`Scénario04_40%`[[20]],
                              comparaison01$`Scénario01_0%`[[26]], comparaison01$`Scénario02_10%`[[26]], comparaison01$`Scénario03_20%`[[26]], comparaison01$`Scénario04_40%`[[26]],
                              comparaison01$`Scénario01_0%`[[32]], comparaison01$`Scénario02_10%`[[32]], comparaison01$`Scénario03_20%`[[32]], comparaison01$`Scénario04_40%`[[32]],
                              comparaison01$`Scénario01_0%`[[38]], comparaison01$`Scénario02_10%`[[38]], comparaison01$`Scénario03_20%`[[38]], comparaison01$`Scénario04_40%`[[38]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 4),
                  Scénario = c(1, 2, 3, 4))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture

dat <- data.frame(x_cov = c(comparaison01$`Scénario01_0%`[[6]], comparaison01$`Scénario02_10%`[[6]], comparaison01$`Scénario03_20%`[[6]], comparaison01$`Scénario04_40%`[[6]],
                            comparaison01$`Scénario01_0%`[[12]], comparaison01$`Scénario02_10%`[[12]], comparaison01$`Scénario03_20%`[[12]], comparaison01$`Scénario04_40%`[[12]],
                            comparaison01$`Scénario01_0%`[[18]], comparaison01$`Scénario02_10%`[[18]], comparaison01$`Scénario03_20%`[[18]], comparaison01$`Scénario04_40%`[[18]],
                            comparaison01$`Scénario01_0%`[[24]], comparaison01$`Scénario02_10%`[[24]], comparaison01$`Scénario03_20%`[[24]], comparaison01$`Scénario04_40%`[[24]],
                            comparaison01$`Scénario01_0%`[[30]], comparaison01$`Scénario02_10%`[[30]], comparaison01$`Scénario03_20%`[[30]], comparaison01$`Scénario04_40%`[[30]],
                            comparaison01$`Scénario01_0%`[[36]], comparaison01$`Scénario02_10%`[[36]], comparaison01$`Scénario03_20%`[[36]], comparaison01$`Scénario04_40%`[[36]],
                            comparaison01$`Scénario01_0%`[[42]], comparaison01$`Scénario02_10%`[[42]], comparaison01$`Scénario03_20%`[[42]], comparaison01$`Scénario04_40%`[[42]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 4),
                  Scénario = c(1, 2, 3, 4))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Biais médian

dat <- data.frame(x_biais = c(comparaison01$`Scénario01_0%`[[1]], comparaison01$`Scénario02_10%`[[1]], comparaison01$`Scénario03_20%`[[1]], comparaison01$`Scénario04_40%`[[1]],
                              comparaison01$`Scénario01_0%`[[7]], comparaison01$`Scénario02_10%`[[7]], comparaison01$`Scénario03_20%`[[7]], comparaison01$`Scénario04_40%`[[7]],
                              comparaison01$`Scénario01_0%`[[13]], comparaison01$`Scénario02_10%`[[13]], comparaison01$`Scénario03_20%`[[13]], comparaison01$`Scénario04_40%`[[13]],
                              comparaison01$`Scénario01_0%`[[19]], comparaison01$`Scénario02_10%`[[19]], comparaison01$`Scénario03_20%`[[19]], comparaison01$`Scénario04_40%`[[19]],
                              comparaison01$`Scénario01_0%`[[25]], comparaison01$`Scénario02_10%`[[25]], comparaison01$`Scénario03_20%`[[25]], comparaison01$`Scénario04_40%`[[25]],
                              comparaison01$`Scénario01_0%`[[31]], comparaison01$`Scénario02_10%`[[31]], comparaison01$`Scénario03_20%`[[31]], comparaison01$`Scénario04_40%`[[31]],
                              comparaison01$`Scénario01_0%`[[37]], comparaison01$`Scénario02_10%`[[37]], comparaison01$`Scénario03_20%`[[37]], comparaison01$`Scénario04_40%`[[37]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 4),
                  Scénario = c(1, 2, 3, 4))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode)) 

################################################################################
################################################################################
## nsim = 5000, Scénarios 05 prévalence de I2  ~ 10%, 15% et 30%, taux de co-infection ~30%-40%

# Tableau de comparaison

comparaison06 <- data.frame(matrix(ncol = 5, 
                                   nrow = 30))

colnames(comparaison06) <- c(
  "Estimation", 
  "Performance",
  "Scénario05_10%",
  "Scénario05_15%",
  "Scénario05_20%"
)

comparaison06$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison06$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison06$`Scénario05_10%` <- c(
  
  0.150884881370104, 0.00261027654528226, 0.0253739136948072, 0.0300021404233505, 0.021,
  0.1438125433141, 0.0025488274676997, 0.0253739136948072, 0.0367713239172938, 0.02,
  0.150505992105703,  0.00281534225959047, 0.0254645805770528, 0.0318342376064884, 0.034,
  0.147404111056842, 0.0025951574942643, 0.0243205342932279, 0.0458889476395196, 0.011,
  0.136551656896362, 0.00964930135723267, 0.0282860070570173, 0.608349891873479, 0.155,
  0.13997310413075, 0.00524103605401349, 0.0248282648979573, 0.312415604406487, 0.076
  
)

comparaison06$`Scénario05_15%` <- c(
  
  0.140429123280203, 0.00249710455610192, 0.0222149461167922, 0.0271294734155219, 0.029,
  0.129221105313994, 0.00243817440796002, 0.0222149461167922, 0.0328697740340048, 0.045,
  0.134370408320079,  0.00267760305174534, 0.0207303320807982, 0.0292850396000047, 0.069,
  0.133380796854627, 0.00251658008156455, 0.0203045004710583, 0.0500133138140319, 0.013,
  0.12659985027029, 0.00770953870299387, 0.0237293512527508, 0.454224430051889, 0.212,
  0.122963364920755, 0.0067282807263848, 0.0218415415582933, 0.556612375414292, 0.131
  
)

comparaison06$`Scénario05_20%` <- c(
  
  0.10502582296353, 0.00229248281275315, 0.0133206138191071, 0.0233585894807925, 0.166,
  0.0895370366610017, 0.00226148026903347, 0.0133206138191071, 0.025357196299331, 0.284,
  0.0900426984752905, 0.0025062913164543, 0.0106114725738499, 0.0273487713593666, 0.339,
  0.125169476540054, 0.0521941435774318, 0.0678093472911656, 1.14083839011641, 0.477,
  0.0947781068980748, 0.00255199206552236, 0.0115323296206397, 0.059023116139601, 0.124,
  0.061031136174096, 0.0156272216457457, 0.019336394006801, 0.846403007930905, 0.492

)
####################### Représentation graphique ###############################
# Biais

dat <- data.frame(x_biais = c(comparaison06$`Scénario05_10%`[[1]], comparaison06$`Scénario05_15%`[[1]], comparaison06$`Scénario05_20%`[[1]],
                               comparaison06$`Scénario05_10%`[[6]], comparaison06$`Scénario05_15%`[[6]], comparaison06$`Scénario05_20%`[[6]],
                               comparaison06$`Scénario05_10%`[[11]], comparaison06$`Scénario05_15%`[[11]], comparaison06$`Scénario05_20%`[[11]],
                               comparaison06$`Scénario05_10%`[[16]], comparaison06$`Scénario05_15%`[[16]], comparaison06$`Scénario05_20%`[[16]],
                               comparaison06$`Scénario05_10%`[[21]], comparaison06$`Scénario05_15%`[[21]], comparaison06$`Scénario05_20%`[[21]],
                               comparaison06$`Scénario05_10%`[[26]], comparaison06$`Scénario05_15%`[[26]], comparaison06$`Scénario05_20%`[[26]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                   Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture

ddat <- data.frame(x_biais = c(comparaison06$`Scénario05_10%`[[5]], comparaison06$`Scénario05_15%`[[5]], comparaison06$`Scénario05_20%`[[5]],
                                comparaison06$`Scénario05_10%`[[10]], comparaison06$`Scénario05_15%`[[10]], comparaison06$`Scénario05_20%`[[10]],
                                comparaison06$`Scénario05_10%`[[15]], comparaison06$`Scénario05_15%`[[15]], comparaison06$`Scénario05_20%`[[15]],
                                comparaison06$`Scénario05_10%`[[20]], comparaison06$`Scénario05_15%`[[20]], comparaison06$`Scénario05_20%`[[20]],
                                comparaison06$`Scénario05_10%`[[25]], comparaison06$`Scénario05_15%`[[25]], comparaison06$`Scénario05_20%`[[25]],
                                comparaison06$`Scénario05_10%`[[30]], comparaison06$`Scénario05_15%`[[30]], comparaison06$`Scénario05_20%`[[30]]),
                    Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                    Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

###############################################################################
## Scénarios avec différentes valeurs de couverture vaccinale et un taux de co-infection entre 30% et 40%

# Tableau de comparaison

comparaison07 <- data.frame(matrix(ncol = 5, 
                                   nrow = 30))

colnames(comparaison07) <- c(
  "Estimation", 
  "Performance",
  "Scénario06_15%",
  "Scénario06_50%",
  "Scénario06_70%"
)

comparaison07$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison07$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison07$`Scénario06_15%` <- c(
  
  0.114785684347962, 0.00577588163677071, 0.018945859086364, 0.0231971843422438, 0.5,
  0.0998165672971633, 0.0058927572201622, 0.018945859086364, 0.0237639608339032, 0.609,
  0.0764219679845992, 0.011162386406034, 0.0169915412102671, 0.224320168725075, 0.815,
  0.0861120675120655, 0.0482415186002367, 0.055608565252839, -1.1192906240961, 0.818,
  0.0371805581801326, 0.0385528610379995, 0.0398967020835477, 0.586869217260683, 0.911,
  -0.00166307065713395, 0.0505673111035422, 0.0505195095964493, -0.77870549573093 , 0.868
  
)

comparaison07$`Scénario06_50%` <- c(
  
  0.113441130946681, 0.00290976236935777, 0.0157757427974504, 0.0201964106881177, 0.22,
  0.0868127367349568, 0.00312807050499393, 0.0157757427974504, 0.0129966430972609, 0.505,
  0.0839960039605492, 0.00338985210140379, 0.010441790930643, 0.0297190116296657, 0.509,
  0.0944630126417745, 0.0100082502027458, 0.0189215027099031, 0.0283567080494619, 0.758,
  0.0726014705401477, 0.0132319989594893, 0.0184897404851217, 0.582380637414296, 0.663,
  0.0570439737517115, 0.0282004660401309, 0.0314262805154767, 0.835088548444357, 0.623
  
)

comparaison07$`Scénario06_70%` <- c(
  
  0.117556895387851, 0.00342870611602003, 0.0172449010631342, 0.0223567732876053, 0.26,
  0.0857554974362157, 0.00382330036249513, 0.0172449010631342, 0.0113772999708081, 0.594,
  0.0773143193308095, 0.0048321098456319, 0.0108047817093726, 0.0397776427849472, 0.647,
  0.0968822730741936, 0.0082804599495293, 0.0176583543256024, 0.0474536977003922, 0.653,
  0.0813216254517906, 0.0101838173883938, 0.0167868403371267, 0.330184033431647, 0.684,
  0.0553564271802355, 0.040736470182651, 0.0437600677426291, 1.50972791862362, 0.627
  
)

View(comparaison07)

####################### Représentation graphique ###############################

# Biais
dat <- data.frame(x_biais = c(comparaison07$`Scénario06_15%`[[1]], comparaison07$`Scénario06_50%`[[1]], comparaison07$`Scénario06_70%`[[1]],
                               comparaison07$`Scénario06_15%`[[6]], comparaison07$`Scénario06_50%`[[6]], comparaison07$`Scénario06_70%`[[6]],
                               comparaison07$`Scénario06_15%`[[11]], comparaison07$`Scénario06_50%`[[11]], comparaison07$`Scénario06_70%`[[11]],
                               comparaison07$`Scénario06_15%`[[16]], comparaison07$`Scénario06_50%`[[16]], comparaison07$`Scénario06_70%`[[16]],
                               comparaison07$`Scénario06_15%`[[21]], comparaison07$`Scénario06_50%`[[21]], comparaison07$`Scénario06_70%`[[21]],
                               comparaison07$`Scénario06_15%`[[26]], comparaison07$`Scénario06_50%`[[26]], comparaison07$`Scénario06_70%`[[26]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                   Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture
dat <- data.frame(x_cov = c(comparaison07$`Scénario06_15%`[[5]], comparaison07$`Scénario06_50%`[[5]], comparaison07$`Scénario06_70%`[[5]],
                             comparaison07$`Scénario06_15%`[[10]], comparaison07$`Scénario06_50%`[[10]], comparaison07$`Scénario06_70%`[[10]],
                             comparaison07$`Scénario06_15%`[[15]], comparaison07$`Scénario06_50%`[[15]], comparaison07$`Scénario06_70%`[[15]],
                             comparaison07$`Scénario06_15%`[[20]], comparaison07$`Scénario06_50%`[[20]], comparaison07$`Scénario06_70%`[[20]],
                             comparaison07$`Scénario06_15%`[[25]], comparaison07$`Scénario06_50%`[[25]], comparaison07$`Scénario06_70%`[[25]],
                             comparaison07$`Scénario06_15%`[[30]], comparaison07$`Scénario06_50%`[[30]], comparaison07$`Scénario06_70%`[[30]]),
                   Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                   Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

###############################################################################
## Scénarios 07 prévalence de I1  ~ 10%, 30% et 50%, taux de co-infection ~30%-40%


# Tableau de comparaison

comparaison08 <- data.frame(matrix(ncol = 5, 
                                   nrow = 30))

colnames(comparaison08) <- c(
  "Estimation", 
  "Performance",
  "Scénario07_10%",
  "Scénario07_20%",
  "Scénario07_50%"
)

comparaison08$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison08$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison08$`Scénario07_10%` <- c(
  
  0.106673624839017, 0.00799997589595763, 0.019371238156357, 0.0389681769237109, 0.603,
  0.084364454614073, 0.00860548401382314, 0.019371238156357, 0.0293893608850529, 0.758,
  0.078908263163963, 0.00928808954300103, 0.0155053154490113, 0.0447567255088245, 0.811,
  -0.0607668379950811, 0.051479124031844, 0.0551202535077326, -1.19530010301766, 0.98,
  0.0333517351021474, 0.0367555631863109, 0.0378311458574484, -0.195292253623931, 0.891,
  -0.0627418319342407, 0.0835603953789433, 0.0874133724580289, 1.45238643308653, 0.851
  
)

comparaison08$`Scénario07_20%` <- c(
  
  0.124329140539258, 0.00183177360896528, 0.0172876770225868, 0.0116810374864103, 0.052,
  0.105189044293054, 0.00181629409710442, 0.0172876770225868, 0.0154678825509317, 0.132,
  0.103542262105581, 0.00191358896704037, 0.0126326754200141, 0.0129435783352396, 0.221,
  0.110591641444496, 0.00189903443264864, 0.014127646555604, 0.0435275832129529, 0.035,
  0.097294580291964, 0.00581728086425926, 0.0152776989375844, 0.38768492267554, 0.346,
  0.0844774715900821, 0.00949269863726164, 0.0166196491448775, 0.836324539254684, 0.325
  
)

comparaison08$`Scénario07_50%` <- c(
  
  0.155235035991219, 0.00176140837796493, 0.0258575633687821, 0.00918210720861067, 0.008,
  0.136648303446157, 0.00173978048134725, 0.0258575633687821, 0.0156709467419737, 0.016,
  0.134853794492636, 0.00191173086759032, 0.0200953650257847, 0.0116063553254991, 0.048,
  0.141210108892723, 0.00179759803375561, 0.0217360952892166, 0.0332508991603645, 0.006,
  0.130695534479945, 0.00568506682487009, 0.0227607044910437, 0.432681825999463, 0.128,
  0.126681809809132, 0.00627518648919158, 0.0223171922392196, 0.639483087298042, 0.08
  
)

View(comparaison08)

####################### Représentation graphique ###############################

# Biais
dat <- data.frame(x_biais = c(comparaison08$`Scénario07_10%`[[1]], comparaison08$`Scénario07_20%`[[1]], comparaison08$`Scénario07_50%`[[1]],
                              comparaison08$`Scénario07_10%`[[6]], comparaison08$`Scénario07_20%`[[6]], comparaison08$`Scénario07_50%`[[6]],
                              comparaison08$`Scénario07_10%`[[11]], comparaison08$`Scénario07_20%`[[11]], comparaison08$`Scénario07_50%`[[11]],
                              comparaison08$`Scénario07_10%`[[16]], comparaison08$`Scénario07_20%`[[16]], comparaison08$`Scénario07_50%`[[16]],
                              comparaison08$`Scénario07_10%`[[21]], comparaison08$`Scénario07_20%`[[21]], comparaison08$`Scénario07_50%`[[21]],
                              comparaison08$`Scénario07_10%`[[26]], comparaison08$`Scénario07_20%`[[26]], comparaison08$`Scénario07_50%`[[26]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture
dat <- data.frame(x_cov = c(comparaison08$`Scénario07_10%`[[5]], comparaison08$`Scénario07_20%`[[5]], comparaison08$`Scénario07_50%`[[5]],
                            comparaison08$`Scénario07_10%`[[10]], comparaison08$`Scénario07_20%`[[10]], comparaison08$`Scénario07_50%`[[10]],
                            comparaison08$`Scénario07_10%`[[15]], comparaison08$`Scénario07_20%`[[15]], comparaison08$`Scénario07_50%`[[15]],
                            comparaison08$`Scénario07_10%`[[20]], comparaison08$`Scénario07_20%`[[20]], comparaison08$`Scénario07_50%`[[20]],
                            comparaison08$`Scénario07_10%`[[25]], comparaison08$`Scénario07_20%`[[25]], comparaison08$`Scénario07_50%`[[25]],
                            comparaison08$`Scénario07_10%`[[30]], comparaison08$`Scénario07_20%`[[30]], comparaison08$`Scénario07_50%`[[30]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))
###############################################################################
## Scénarios 09 prévalence de W1  ~ 10%, 15% et 20%, taux de co-infection ~30%-40%


# Tableau de comparaison

comparaison09 <- data.frame(matrix(ncol = 5, 
                                   nrow = 30))

colnames(comparaison09) <- c(
  "Estimation", 
  "Performance",
  "Scénario09_10%",
  "Scénario09_15%",
  "Scénario09_20%"
)

comparaison09$Estimation <- c("RegLog", "-", "-", "-", "-", "IPW", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-")

comparaison09$Performance <- rep(c("Biais", "Variance", "MSE", "Précision", "%Cov"), 6)

comparaison09$`Scénario09_10%` <- c(
  
  0.165179660374365, 0.00282576029786419, 0.0301072547389569, 0.0177207992444735, 0.016,
  0.142443047454319, 0.00268513584371904, 0.0301072547389569, 0.0167132092487773, 0.066,
  0.13912192754701, 0.00286631584992144, 0.0222183602584669, 0.0167844280214891, 0.121,
  0.153126261887726, 0.00979862380768406, 0.0332364772635847, 0.0404444085787871, 0.465,
  0.128972639384344, 0.0119658325648453, 0.0285878084420444, 0.563806166034017, 0.283,
  -0.0841313562345644, 0.0294098984537943, 0.0364585736572077, 1.3487946440389, 0.355
  
)

comparaison09$`Scénario09_15%` <- c(
  
  0.196272852751467, 0.00235204774254418, 0.0408727284220008, 0.0119191357740873, 0.001,
  0.176015118415779, 0.00235412468397975, 0.0408727284220008, 0.0162077255678053, 0.006,
  0.17307294777928, 0.00255019007953944, 0.0325018851424692, 0.0132918507089151, 0.017,
  0.182886261414897, 0.00265804978344263, 0.036102776347977, 0.0388477114718815, 0.048,
  0.164152124944987, 0.00779559015491422, 0.0347337146887138, 0.350179353478576, 0.116,
  0.150737755370414, 0.018449657641898, 0.041153078878367, 0.864217254207267, 0.134 
  
)

comparaison09$`Scénario09_20%` <- c(
  
  0.234877793407975, 0.00253954395291001, 0.0577045822451564, 0.0134101440356874, 0,
  0.216378511934917, 0.00248830623115789, 0.0577045822451564, 0.0179464566208088, 0,
  0.21244977417758, 0.00260385852604288, 0.0477361612156218, 0.0129443619646356, 0,
  0.22139179082562, 0.00258998712974172, 0.0516017221875871, 0.0417855853198331, 0.004,
  0.206238900631201, 0.00755724818612579, 0.050084175071506, 0.458534613959005, 0.058,
  0.193742774410024, 0.0131676223579444, 0.0506907173716801, 0.741852846879047, 0.075
  
)

View(comparaison09)

####################### Représentation graphique ###############################

# Biais
dat <- data.frame(x_biais = c(comparaison09$`Scénario09_10%`[[1]], comparaison09$`Scénario09_15%`[[1]], comparaison09$`Scénario09_20%`[[1]],
                              comparaison09$`Scénario09_10%`[[6]], comparaison09$`Scénario09_15%`[[6]], comparaison09$`Scénario09_20%`[[6]],
                              comparaison09$`Scénario09_10%`[[11]], comparaison09$`Scénario09_15%`[[11]], comparaison09$`Scénario09_20%`[[11]],
                              comparaison09$`Scénario09_10%`[[16]], comparaison09$`Scénario09_15%`[[16]], comparaison09$`Scénario09_20%`[[16]],
                              comparaison09$`Scénario09_10%`[[21]], comparaison09$`Scénario09_15%`[[21]], comparaison09$`Scénario09_20%`[[21]],
                              comparaison09$`Scénario09_10%`[[26]], comparaison09$`Scénario09_15%`[[26]], comparaison09$`Scénario09_20%`[[26]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))

# Couverture
dat <- data.frame(x_cov = c(comparaison09$`Scénario09_10%`[[5]], comparaison09$`Scénario09_15%`[[5]], comparaison09$`Scénario09_20%`[[5]],
                            comparaison09$`Scénario09_10%`[[10]], comparaison09$`Scénario09_15%`[[10]], comparaison09$`Scénario09_20%`[[10]],
                            comparaison09$`Scénario09_10%`[[15]], comparaison09$`Scénario09_15%`[[15]], comparaison09$`Scénario09_20%`[[15]],
                            comparaison09$`Scénario09_10%`[[20]], comparaison09$`Scénario09_15%`[[20]], comparaison09$`Scénario09_20%`[[20]],
                            comparaison09$`Scénario09_10%`[[25]], comparaison09$`Scénario09_15%`[[25]], comparaison09$`Scénario09_20%`[[25]],
                            comparaison09$`Scénario09_10%`[[30]], comparaison09$`Scénario09_15%`[[30]], comparaison09$`Scénario09_20%`[[30]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat) +
  aes(x = Scénario, y = x_cov) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))