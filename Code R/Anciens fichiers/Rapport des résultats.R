# Chargement des librairies nécessaires

library(ggplot2)
library(patchwork)
library(kableExtra)
library("gridExtra")
library("cowplot")
library(ggpubr)
library(Rmisc)

################################################################################
################################################################################
## Scénarios avec différents taux de co-infection : 0%, 10%, 20% et 40%
# nsim = 1000, Scénarios 01, 02, 03 et 04

# Chargement des résultats bruts

## Taille d'échantillo = 1000 ##

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario01.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario02.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario03.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario04.RData")

## Taille d'échantillo = 5000 ##

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario01_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario02_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario03_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario04_5000.RData")

# Tableau de comparaison 

comparaison01 <- data.frame(matrix(ncol = 10, 
                                   nrow = 42))

colnames(comparaison01) <- c(
  "Estimation", 
  "Performance",
  
  "Scénario01_0%",
  "Scénario01_0%_5000",
  
  "Scénario02_10%",
  "Scénario02_10%_5000",
  
  "Scénario03_20%",
  "Scénario03_20%_5000",
  
  "Scénario04_40%",
  "Scénario04_40%_5000"
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

comparaison01$`Scénario01_0%_5000` <- c(
  
  (median(resultats_1_5000$RRc) - mean(l_vraiRRc_1_5000)), Tab01_1_5000$Autres, Tab01_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_1_5000$RRm) - mean(l_vraiRRm_1_5000)), Tab02_1_5000$Autres, Tab02_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_1_5000$RRm_RF) - mean(l_vraiRRm_1_5000)), Tab03_1_5000$Autres[1:4], Tab03_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_1_5000$RRm_Lasso) - mean(l_vraiRRm_1_5000)),Tab03_1_5000$Autres[5:8], Tab03_1_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_1_5000$RRm_Mars) - mean(l_vraiRRm_1_5000)), Tab03_1_5000$Autres[9:12], Tab03_1_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_1_5000$RRm_RN) - mean(l_vraiRRm_1_5000)), Tab03_1_5000$Autres[13:16], Tab03_1_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_1_5000$RRm_GLM) - mean(l_vraiRRm_1_5000)), Tab03_1_5000$Autres[17:20], Tab03_1_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario02_10%_5000` <- c(
  
  (median(resultats_2_5000$RRc) - mean(l_vraiRRc_2_5000)), Tab01_2_5000$Autres, Tab01_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_2_5000$RRm) - mean(l_vraiRRm_2_5000)), Tab02_2_5000$Autres, Tab02_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_2_5000$RRm_RF) - mean(l_vraiRRm_2_5000)), Tab03_2_5000$Autres[1:4], Tab03_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_2_5000$RRm_Lasso) - mean(l_vraiRRm_2_5000)),Tab03_2_5000$Autres[5:8], Tab03_2_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_2_5000$RRm_Mars) - mean(l_vraiRRm_2_5000)), Tab03_2_5000$Autres[9:12], Tab03_2_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_2_5000$RRm_RN) - mean(l_vraiRRm_2_5000)), Tab03_2_5000$Autres[13:16], Tab03_2_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_2_5000$RRm_GLM) - mean(l_vraiRRm_2_5000)), Tab03_2_5000$Autres[17:20], Tab03_2_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario03_20%_5000` <- c(
  
  (median(resultats_3_5000$RRc) - mean(l_vraiRRc_3_5000)), Tab01_3_5000$Autres, Tab01_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_3_5000$RRm) - mean(l_vraiRRm_3_5000)), Tab02_3_5000$Autres, Tab02_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_3_5000$RRm_RF) - mean(l_vraiRRm_3_5000)), Tab03_3_5000$Autres[1:4], Tab03_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_3_5000$RRm_Lasso) - mean(l_vraiRRm_3_5000)),Tab03_3_5000$Autres[5:8], Tab03_3_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_3_5000$RRm_Mars) - mean(l_vraiRRm_3_5000)), Tab03_3_5000$Autres[9:12], Tab03_3_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_3_5000$RRm_RN) - mean(l_vraiRRm_3_5000)), Tab03_3_5000$Autres[13:16], Tab03_3_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_3_5000$RRm_GLM) - mean(l_vraiRRm_3_5000)), Tab03_3_5000$Autres[17:20], Tab03_3_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison01$`Scénario04_40%_5000` <- c(
  
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

comparaison01$`Scénario01_0%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario01_0%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison01$`Scénario02_10%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario02_10%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison01$`Scénario03_20%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario03_20%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison01$`Scénario04_40%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison01$`Scénario04_40%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])


View(comparaison01)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison01$`Scénario01_0%`[[2]], comparaison01$`Scénario02_10%`[[2]], comparaison01$`Scénario03_20%`[[2]], comparaison01$`Scénario04_40%`[[2]],
                              comparaison01$`Scénario01_0%_5000`[[2]], comparaison01$`Scénario02_10%_5000`[[2]], comparaison01$`Scénario03_20%_5000`[[2]], comparaison01$`Scénario04_40%_5000`[[2]],
                              
                              comparaison01$`Scénario01_0%`[[8]], comparaison01$`Scénario02_10%`[[8]], comparaison01$`Scénario03_20%`[[8]], comparaison01$`Scénario04_40%`[[8]],
                              comparaison01$`Scénario01_0%_5000`[[8]], comparaison01$`Scénario02_10%_5000`[[8]], comparaison01$`Scénario03_20%_5000`[[8]], comparaison01$`Scénario04_40%_5000`[[8]],
                              
                              comparaison01$`Scénario01_0%`[[14]], comparaison01$`Scénario02_10%`[[14]], comparaison01$`Scénario03_20%`[[14]], comparaison01$`Scénario04_40%`[[14]],
                              comparaison01$`Scénario01_0%_5000`[[14]], comparaison01$`Scénario02_10%_5000`[[14]], comparaison01$`Scénario03_20%_5000`[[14]], comparaison01$`Scénario04_40%_5000`[[14]],
                              
                              comparaison01$`Scénario01_0%`[[20]], comparaison01$`Scénario02_10%`[[20]], comparaison01$`Scénario03_20%`[[20]], comparaison01$`Scénario04_40%`[[20]],
                              comparaison01$`Scénario01_0%_5000`[[20]], comparaison01$`Scénario02_10%_5000`[[20]], comparaison01$`Scénario03_20%_5000`[[20]], comparaison01$`Scénario04_40%_5000`[[20]],
                              
                              comparaison01$`Scénario01_0%`[[26]], comparaison01$`Scénario02_10%`[[26]], comparaison01$`Scénario03_20%`[[26]], comparaison01$`Scénario04_40%`[[26]],
                              comparaison01$`Scénario01_0%_5000`[[26]], comparaison01$`Scénario02_10%_5000`[[26]], comparaison01$`Scénario03_20%_5000`[[26]], comparaison01$`Scénario04_40%_5000`[[26]],
                              
                              comparaison01$`Scénario01_0%`[[32]], comparaison01$`Scénario02_10%`[[32]], comparaison01$`Scénario03_20%`[[32]], comparaison01$`Scénario04_40%`[[32]],
                              comparaison01$`Scénario01_0%_5000`[[32]], comparaison01$`Scénario02_10%_5000`[[32]], comparaison01$`Scénario03_20%_5000`[[32]], comparaison01$`Scénario04_40%_5000`[[32]],
                              
                              comparaison01$`Scénario01_0%`[[38]], comparaison01$`Scénario02_10%`[[38]], comparaison01$`Scénario03_20%`[[38]], comparaison01$`Scénario04_40%`[[38]],
                              comparaison01$`Scénario01_0%_5000`[[38]], comparaison01$`Scénario02_10%_5000`[[38]], comparaison01$`Scénario03_20%_5000`[[38]], comparaison01$`Scénario04_40%_5000`[[38]]),
                  
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 8),
                  Taille = rep(c("1000", "5000"), each = 4),
                  Scénario = c(1, 2, 3, 4))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)
dat$Taille <- as.factor(dat$Taille)


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                      group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  

  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)
  
  theme_bw() 

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


ggplot(dat, aes(x = Scénario, y = x_cov, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

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


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

################################################################################
################################################################################
# Scénarios avec différentes valeurs de prévalence de I2 et un taux de co-infection entre 30% et 40%
# nsim = 1000, prévalence de I2  ~ 10%, 15% et 30%

## Chargement des résultats bruts

## Taille d'échantillo = 1000 ##

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario05_10%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario05_15%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario05_30%.RData")

## Taille d'échantillo = 5000 ##

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario05_10%_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario05_15%_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario05_30%_5000.RData")

comparaison02 <- data.frame(matrix(ncol = 8, 
                                   nrow = 42))

colnames(comparaison02) <- c(
  "Estimation", 
  "Performance",
  "Scénario05_10%",
  "Scénario05_15%",
  "Scénario05_30%",
  "Scénario05_10%_5000",
  "Scénario05_15%_5000",
  "Scénario05_30%_5000"
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

comparaison02$`Scénario05_10%_5000` <- c(
  
  (median(resultats_5_1_5000$RRc) - mean(l_vraiRRc_5_1_5000)), Tab01_5_1_5000$Autres, Tab01_5_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_5_1_5000$RRm) - mean(l_vraiRRm_5_1_5000)), Tab02_5_1_5000$Autres, Tab02_5_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_1_5000$RRm_RF) - mean(l_vraiRRm_5_1_5000)), Tab03_5_1_5000$Autres[1:4], Tab03_5_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_1_5000$RRm_Lasso) - mean(l_vraiRRm_5_1_5000)),Tab03_5_1_5000$Autres[5:8], Tab03_5_1_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_5_1_5000$RRm_Mars) - mean(l_vraiRRm_5_1_5000)), Tab03_5_1_5000$Autres[9:12], Tab03_5_1_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_5_1_5000$RRm_RN) - mean(l_vraiRRm_5_1_5000)), Tab03_5_1_5000$Autres[13:16], Tab03_5_1_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_5_1_5000$RRm_GLM) - mean(l_vraiRRm_5_1_5000)), Tab03_5_1_5000$Autres[17:20], Tab03_5_1_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison02$`Scénario05_15%_5000` <- c(
  
  (median(resultats_5_2_5000$RRc) - mean(l_vraiRRc_5_2_5000)), Tab01_5_2_5000$Autres, Tab01_5_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_5_2_5000$RRm) - mean(l_vraiRRm_5_2_5000)), Tab02_5_2_5000$Autres, Tab02_5_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_2_5000$RRm_RF) - mean(l_vraiRRm_5_2_5000)), Tab03_5_2_5000$Autres[1:4], Tab03_5_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_2_5000$RRm_Lasso) - mean(l_vraiRRm_5_2_5000)),Tab03_5_2_5000$Autres[5:8], Tab03_5_2_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_5_2_5000$RRm_Mars) - mean(l_vraiRRm_5_2_5000)), Tab03_5_2_5000$Autres[9:12], Tab03_5_2_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_5_2_5000$RRm_RN) - mean(l_vraiRRm_5_2_5000)), Tab03_5_2_5000$Autres[13:16], Tab03_5_2_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_5_2_5000$RRm_GLM) - mean(l_vraiRRm_5_2_5000)), Tab03_5_2_5000$Autres[17:20], Tab03_5_2_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison02$`Scénario05_30%_5000` <- c(
  
  (median(resultats_5_3_5000$RRc) - mean(l_vraiRRc_5_3_5000)), Tab01_5_3_5000$Autres, Tab01_5_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_5_3_5000$RRm) - mean(l_vraiRRm_5_3_5000)), Tab02_5_3_5000$Autres, Tab02_5_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_3_5000$RRm_RF) - mean(l_vraiRRm_5_3_5000)), Tab03_5_3_5000$Autres[1:4], Tab03_5_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_5_3_5000$RRm_Lasso) - mean(l_vraiRRm_5_3_5000)),Tab03_5_3_5000$Autres[5:8], Tab03_5_3_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_5_3_5000$RRm_Mars) - mean(l_vraiRRm_5_3_5000)), Tab03_5_3_5000$Autres[9:12], Tab03_5_3_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_5_3_5000$RRm_RN) - mean(l_vraiRRm_5_3_5000)), Tab03_5_3_5000$Autres[13:16], Tab03_5_3_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_5_3_5000$RRm_GLM) - mean(l_vraiRRm_5_3_5000)), Tab03_5_3_5000$Autres[17:20], Tab03_5_3_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison02$`Scénario05_10%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison02$`Scénario05_10%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison02$`Scénario05_15%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison02$`Scénario05_15%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison02$`Scénario05_30%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison02$`Scénario05_30%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

comparaison02$`Scénario05_10%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison02$`Scénario05_10%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison02$`Scénario05_15%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison02$`Scénario05_15%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison02$`Scénario05_30%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison02$`Scénario05_30%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison02)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison02$`Scénario05_10%`[[2]], comparaison02$`Scénario05_15%`[[2]], comparaison02$`Scénario05_30%`[[2]],
                              comparaison02$`Scénario05_10%_5000`[[2]], comparaison02$`Scénario05_15%_5000`[[2]], comparaison02$`Scénario05_30%_5000`[[2]],
                              
                              comparaison02$`Scénario05_10%`[[8]], comparaison02$`Scénario05_15%`[[8]], comparaison02$`Scénario05_30%`[[8]],
                              comparaison02$`Scénario05_10%_5000`[[8]], comparaison02$`Scénario05_15%_5000`[[8]], comparaison02$`Scénario05_30%_5000`[[8]],
                              
                              comparaison02$`Scénario05_10%`[[14]], comparaison02$`Scénario05_15%`[[14]], comparaison02$`Scénario05_30%`[[14]],
                              comparaison02$`Scénario05_10%_5000`[[14]], comparaison02$`Scénario05_15%_5000`[[14]], comparaison02$`Scénario05_30%_5000`[[14]],
                              
                              comparaison02$`Scénario05_10%`[[20]], comparaison02$`Scénario05_15%`[[20]], comparaison02$`Scénario05_30%`[[20]],
                              comparaison02$`Scénario05_10%_5000`[[20]], comparaison02$`Scénario05_15%_5000`[[20]], comparaison02$`Scénario05_30%_5000`[[20]],
                              
                              comparaison02$`Scénario05_10%`[[26]], comparaison02$`Scénario05_15%`[[26]], comparaison02$`Scénario05_30%`[[26]],
                              comparaison02$`Scénario05_10%_5000`[[26]], comparaison02$`Scénario05_15%_5000`[[26]], comparaison02$`Scénario05_30%_5000`[[26]],
                              
                              comparaison02$`Scénario05_10%`[[32]], comparaison02$`Scénario05_15%`[[32]], comparaison02$`Scénario05_30%`[[32]],
                              comparaison02$`Scénario05_10%_5000`[[32]], comparaison02$`Scénario05_15%_5000`[[32]], comparaison02$`Scénario05_30%_5000`[[32]],
                              
                              comparaison02$`Scénario05_10%`[[38]], comparaison02$`Scénario05_15%`[[38]], comparaison02$`Scénario05_30%`[[38]],
                              comparaison02$`Scénario05_10%_5000`[[38]], comparaison02$`Scénario05_15%_5000`[[38]], comparaison02$`Scénario05_30%_5000`[[38]]),

                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 6),
                  Taille = rep(c("1000", "5000"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

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


ggplot(dat, aes(x = Scénario, y = x_cov, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

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


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

################################################################################
################################################################################
# Scénarios avec différentes valeurs de couverture vaccinale et un taux de co-infection entre 30% et 40%
# nsim = 1000, Scénarios 06 couverture vaccinale  ~ 50%, 70% et 85%

## Chargement des résultats bruts

## Taille d'échantillon n = 1000 ##

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario06_50%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario06_70%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario06_85%.RData")

## Taille d'échantillon n  = 5000 ##

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario06_50%_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario06_70%_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario06_85%_5000.RData")

comparaison03 <- data.frame(matrix(ncol = 8, 
                                   nrow = 42))

colnames(comparaison03) <- c(
  "Estimation", 
  "Performance",
  "Scénario06_50%",
  "Scénario06_70%",
  "Scénario06_85%",
  "Scénario06_50%_5000",
  "Scénario06_70%_5000",
  "Scénario06_85%_5000"
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

comparaison03$`Scénario06_50%_5000` <- c(
  
  (median(resultats_6_1_5000$RRc) - mean(l_vraiRRc_6_1_5000)), Tab01_6_1_5000$Autres, Tab01_6_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_6_1_5000$RRm) - mean(l_vraiRRm_6_1_5000)), Tab02_6_1_5000$Autres, Tab02_6_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_1_5000$RRm_RF) - mean(l_vraiRRm_6_1_5000)), Tab03_6_1_5000$Autres[1:4], Tab03_6_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_1_5000$RRm_Lasso) - mean(l_vraiRRm_6_1_5000)),Tab03_6_1_5000$Autres[5:8], Tab03_6_1_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_6_1_5000$RRm_Mars) - mean(l_vraiRRm_6_1_5000)), Tab03_6_1_5000$Autres[9:12], Tab03_6_1_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_6_1_5000$RRm_RN) - mean(l_vraiRRm_6_1_5000)), Tab03_6_1_5000$Autres[13:16], Tab03_6_1_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_6_1_5000$RRm_GLM) - mean(l_vraiRRm_6_1_5000)), Tab03_6_1_5000$Autres[17:20], Tab03_6_1_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison03$`Scénario06_70%_5000` <- c(
  
  (median(resultats_6_2_5000$RRc) - mean(l_vraiRRc_6_2_5000)), Tab01_6_2_5000$Autres, Tab01_6_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_6_2_5000$RRm) - mean(l_vraiRRm_6_2_5000)), Tab02_6_2_5000$Autres, Tab02_6_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_2_5000$RRm_RF) - mean(l_vraiRRm_6_2_5000)), Tab03_6_2_5000$Autres[1:4], Tab03_6_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_2_5000$RRm_Lasso) - mean(l_vraiRRm_6_2_5000)),Tab03_6_2_5000$Autres[5:8], Tab03_6_2_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_6_2_5000$RRm_Mars) - mean(l_vraiRRm_6_2_5000)), Tab03_6_2_5000$Autres[9:12], Tab03_6_2_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_6_2_5000$RRm_RN) - mean(l_vraiRRm_6_2_5000)), Tab03_6_2_5000$Autres[13:16], Tab03_6_2_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_6_2_5000$RRm_GLM) - mean(l_vraiRRm_6_2_5000)), Tab03_6_2_5000$Autres[17:20], Tab03_6_2_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison03$`Scénario06_85%_5000` <- c(
  
  (median(resultats_6_3_5000$RRc) - mean(l_vraiRRc_6_3_5000)), Tab01_6_3_5000$Autres, Tab01_6_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_6_3_5000$RRm) - mean(l_vraiRRm_6_3_5000)), Tab02_6_3_5000$Autres, Tab02_6_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_3_5000$RRm_RF) - mean(l_vraiRRm_6_3_5000)), Tab03_6_3_5000$Autres[1:4], Tab03_6_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_6_3_5000$RRm_Lasso) - mean(l_vraiRRm_6_3_5000)),Tab03_6_3_5000$Autres[5:8], Tab03_6_3_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_6_3_5000$RRm_Mars) - mean(l_vraiRRm_6_3_5000)), Tab03_6_3_5000$Autres[9:12], Tab03_6_3_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_6_3_5000$RRm_RN) - mean(l_vraiRRm_6_3_5000)), Tab03_6_3_5000$Autres[13:16], Tab03_6_3_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_6_3_5000$RRm_GLM) - mean(l_vraiRRm_6_3_5000)), Tab03_6_3_5000$Autres[17:20], Tab03_6_3_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison03$`Scénario06_50%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison03$`Scénario06_50%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison03$`Scénario06_70%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison03$`Scénario06_70%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison03$`Scénario06_85%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison03$`Scénario06_85%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

comparaison03$`Scénario06_50%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison03$`Scénario06_50%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison03$`Scénario06_70%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison03$`Scénario06_70%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison03$`Scénario06_85%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison03$`Scénario06_85%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison03)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison03$`Scénario06_50%`[[2]], comparaison03$`Scénario06_70%`[[2]], comparaison03$`Scénario06_85%`[[2]],
                              comparaison03$`Scénario06_50%_5000`[[2]], comparaison03$`Scénario06_70%_5000`[[2]], comparaison03$`Scénario06_85%_5000`[[2]],
                              
                              comparaison03$`Scénario06_50%`[[8]], comparaison03$`Scénario06_70%`[[8]], comparaison03$`Scénario06_85%`[[8]],
                              comparaison03$`Scénario06_50%_5000`[[8]], comparaison03$`Scénario06_70%_5000`[[8]], comparaison03$`Scénario06_85%_5000`[[8]],
                              
                              comparaison03$`Scénario06_50%`[[14]], comparaison03$`Scénario06_70%`[[14]], comparaison03$`Scénario06_85%`[[14]],
                              comparaison03$`Scénario06_50%_5000`[[14]], comparaison03$`Scénario06_70%_5000`[[14]], comparaison03$`Scénario06_85%_5000`[[14]],
                              
                              comparaison03$`Scénario06_50%`[[20]], comparaison03$`Scénario06_70%`[[20]], comparaison03$`Scénario06_85%`[[20]],
                              comparaison03$`Scénario06_50%_5000`[[20]], comparaison03$`Scénario06_70%_5000`[[20]], comparaison03$`Scénario06_85%_5000`[[20]],
                              
                              comparaison03$`Scénario06_50%`[[26]], comparaison03$`Scénario06_70%`[[26]], comparaison03$`Scénario06_85%`[[26]],
                              comparaison03$`Scénario06_50%_5000`[[26]], comparaison03$`Scénario06_70%_5000`[[26]], comparaison03$`Scénario06_85%_5000`[[26]],
                              
                              comparaison03$`Scénario06_50%`[[32]], comparaison03$`Scénario06_70%`[[32]], comparaison03$`Scénario06_85%`[[32]],
                              comparaison03$`Scénario06_50%_5000`[[32]], comparaison03$`Scénario06_70%_5000`[[32]], comparaison03$`Scénario06_85%_5000`[[32]],
                              
                              comparaison03$`Scénario06_50%`[[38]], comparaison03$`Scénario06_70%`[[38]], comparaison03$`Scénario06_85%`[[38]],
                              comparaison03$`Scénario06_50%_5000`[[38]], comparaison03$`Scénario06_70%_5000`[[38]], comparaison03$`Scénario06_85%_5000`[[38]]),
                  
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 6),
                  Taille = rep(c("1000", "5000"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

# Couverture

dat <- data.frame(x_cov = c(comparaison03$`Scénario06_50%`[[6]], comparaison03$`Scénario06_70%`[[6]], comparaison03$`Scénario06_85%`[[6]],
                            comparaison03$`Scénario06_50%`[[12]], comparaison03$`Scénario06_70%`[[12]], comparaison03$`Scénario06_85%`[[12]],
                            comparaison03$`Scénario06_50%`[[18]], comparaison03$`Scénario06_70%`[[18]], comparaison03$`Scénario06_85%`[[18]],
                            comparaison03$`Scénario06_50%`[[24]], comparaison03$`Scénario06_70%`[[24]], comparaison03$`Scénario06_85%`[[24]],
                            comparaison03$`Scénario06_50%`[[30]], comparaison03$`Scénario06_70%`[[30]], comparaison03$`Scénario06_85%`[[30]],
                            comparaison03$`Scénario06_50%`[[36]], comparaison03$`Scénario06_70%`[[36]], comparaison03$`Scénario06_85%`[[36]],
                            comparaison03$`Scénario06_50%`[[42]], comparaison03$`Scénario06_70%`[[42]], comparaison03$`Scénario06_85%`[[42]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_cov, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

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


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

################################################################################
################################################################################
# Scénarios avec différentes valeurs de prévalence de I1 et un taux de co-infection entre 30% et 40%
# nsim = 1000, Scénarios 07 prévalence de I1  ~ 10%, 30% et 50%

## Chargement des résultats bruts

## Taille d'échantillon n  = 1000 ##

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario07_10%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario07_30%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario07_50%.RData")

## Taille d'échantillon n  = 5000 ##

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario07_10%_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario07_30%_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario07_50%_5000.RData")


comparaison04 <- data.frame(matrix(ncol = 8, 
                                   nrow = 42))

colnames(comparaison04) <- c(
  "Estimation", 
  "Performance",
  "Scénario07_10%",
  "Scénario07_20%",
  "Scénario07_50%",
  "Scénario07_10%_5000",
  "Scénario07_20%_5000",
  "Scénario07_50%_5000"
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

comparaison04$`Scénario07_10%_5000` <- c(
  
  (median(resultats_7_1_5000$RRc) - mean(l_vraiRRc_7_1_5000)), Tab01_7_1_5000$Autres, Tab01_7_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_7_1_5000$RRm) - mean(l_vraiRRm_7_1_5000)), Tab02_7_1_5000$Autres, Tab02_7_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_1_5000$RRm_RF) - mean(l_vraiRRm_7_1_5000)), Tab03_7_1_5000$Autres[1:4], Tab03_7_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_1_5000$RRm_Lasso) - mean(l_vraiRRm_7_1_5000)),Tab03_7_1_5000$Autres[5:8], Tab03_7_1_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_7_1_5000$RRm_Mars) - mean(l_vraiRRm_7_1_5000)), Tab03_7_1_5000$Autres[9:12], Tab03_7_1_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_7_1_5000$RRm_RN) - mean(l_vraiRRm_7_1_5000)), Tab03_7_1_5000$Autres[13:16], Tab03_7_1_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_7_1_5000$RRm_GLM) - mean(l_vraiRRm_7_1_5000)), Tab03_7_1_5000$Autres[17:20], Tab03_7_1_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison04$`Scénario07_20%_5000` <- c(
  
  (median(resultats_7_2_5000$RRc) - mean(l_vraiRRc_7_2_5000)), Tab01_7_2_5000$Autres, Tab01_7_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_7_2_5000$RRm) - mean(l_vraiRRm_7_2_5000)), Tab02_7_2_5000$Autres, Tab02_7_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_2_5000$RRm_RF) - mean(l_vraiRRm_7_2_5000)), Tab03_7_2_5000$Autres[1:4], Tab03_7_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_2_5000$RRm_Lasso) - mean(l_vraiRRm_7_2_5000)),Tab03_7_2_5000$Autres[5:8], Tab03_7_2_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_7_2_5000$RRm_Mars) - mean(l_vraiRRm_7_2_5000)), Tab03_7_2_5000$Autres[9:12], Tab03_7_2_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_7_2_5000$RRm_RN) - mean(l_vraiRRm_7_2_5000)), Tab03_7_2_5000$Autres[13:16], Tab03_7_2_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_7_2_5000$RRm_GLM) - mean(l_vraiRRm_7_2_5000)), Tab03_7_2_5000$Autres[17:20], Tab03_7_2_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison04$`Scénario07_50%_5000` <- c(
  
  (median(resultats_7_3_5000$RRc) - mean(l_vraiRRc_7_3_5000)), Tab01_7_3_5000$Autres, Tab01_7_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_7_3_5000$RRm) - mean(l_vraiRRm_7_3_5000)), Tab02_7_3_5000$Autres, Tab02_7_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_3_5000$RRm_RF) - mean(l_vraiRRm_7_3_5000)), Tab03_7_3_5000$Autres[1:4], Tab03_7_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_7_3_5000$RRm_Lasso) - mean(l_vraiRRm_7_3_5000)),Tab03_7_3_5000$Autres[5:8], Tab03_7_3_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_7_3_5000$RRm_Mars) - mean(l_vraiRRm_7_3_5000)), Tab03_7_3_5000$Autres[9:12], Tab03_7_3_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_7_3_5000$RRm_RN) - mean(l_vraiRRm_7_3_5000)), Tab03_7_3_5000$Autres[13:16], Tab03_7_3_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_7_3_5000$RRm_GLM) - mean(l_vraiRRm_7_3_5000)), Tab03_7_3_5000$Autres[17:20], Tab03_7_3_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison04$`Scénario07_10%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison04$`Scénario07_10%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison04$`Scénario07_20%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison04$`Scénario07_20%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison04$`Scénario07_50%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison04$`Scénario07_50%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

comparaison04$`Scénario07_10%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison04$`Scénario07_10%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison04$`Scénario07_20%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison04$`Scénario07_20%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison04$`Scénario07_50%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison04$`Scénario07_50%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison04)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison04$`Scénario07_10%`[[2]], comparaison04$`Scénario07_50%`[[2]], comparaison04$`Scénario07_50%`[[2]],
                              comparaison04$`Scénario07_10%_5000`[[2]], comparaison04$`Scénario07_50%_5000`[[2]], comparaison04$`Scénario07_50%_5000`[[2]],
                              
                              comparaison04$`Scénario07_10%`[[8]], comparaison04$`Scénario07_50%`[[8]], comparaison04$`Scénario07_50%`[[8]],
                              comparaison04$`Scénario07_10%_5000`[[8]], comparaison04$`Scénario07_50%_5000`[[8]], comparaison04$`Scénario07_50%_5000`[[8]],
                              
                              comparaison04$`Scénario07_10%`[[14]], comparaison04$`Scénario07_50%`[[14]], comparaison04$`Scénario07_50%`[[14]],
                              comparaison04$`Scénario07_10%_5000`[[14]], comparaison04$`Scénario07_50%_5000`[[14]], comparaison04$`Scénario07_50%_5000`[[14]],
                              
                              comparaison04$`Scénario07_10%`[[20]], comparaison04$`Scénario07_50%`[[20]], comparaison04$`Scénario07_50%`[[20]],
                              comparaison04$`Scénario07_10%_5000`[[20]], comparaison04$`Scénario07_50%_5000`[[20]], comparaison04$`Scénario07_50%_5000`[[20]],
                              
                              comparaison04$`Scénario07_10%`[[26]], comparaison04$`Scénario07_50%`[[26]], comparaison04$`Scénario07_50%`[[26]],
                              comparaison04$`Scénario07_10%_5000`[[26]], comparaison04$`Scénario07_50%_5000`[[26]], comparaison04$`Scénario07_50%_5000`[[26]],
                              
                              comparaison04$`Scénario07_10%`[[32]], comparaison04$`Scénario07_50%`[[32]], comparaison04$`Scénario07_50%`[[32]],
                              comparaison04$`Scénario07_10%_5000`[[32]], comparaison04$`Scénario07_50%_5000`[[32]], comparaison04$`Scénario07_50%_5000`[[32]],
                              
                              comparaison04$`Scénario07_10%`[[38]], comparaison04$`Scénario07_50%`[[38]], comparaison04$`Scénario07_50%`[[38]],
                              comparaison04$`Scénario07_10%_5000`[[38]], comparaison04$`Scénario07_50%_5000`[[38]], comparaison04$`Scénario07_50%_5000`[[38]]),
                  
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 6),
                  Taille = rep(c("1000", "5000"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

# Couverture

dat <- data.frame(x_cov = c(comparaison04$`Scénario07_10%`[[6]], comparaison04$`Scénario07_50%`[[6]], comparaison04$`Scénario07_50%`[[6]],
                            comparaison04$`Scénario07_10%`[[12]], comparaison04$`Scénario07_50%`[[12]], comparaison04$`Scénario07_50%`[[12]],
                            comparaison04$`Scénario07_10%`[[18]], comparaison04$`Scénario07_50%`[[18]], comparaison04$`Scénario07_50%`[[18]],
                            comparaison04$`Scénario07_10%`[[24]], comparaison04$`Scénario07_50%`[[24]], comparaison04$`Scénario07_50%`[[24]],
                            comparaison04$`Scénario07_10%`[[30]], comparaison04$`Scénario07_50%`[[30]], comparaison04$`Scénario07_50%`[[30]],
                            comparaison04$`Scénario07_10%`[[36]], comparaison04$`Scénario07_50%`[[36]], comparaison04$`Scénario07_50%`[[36]],
                            comparaison04$`Scénario07_10%`[[42]], comparaison04$`Scénario07_50%`[[42]], comparaison04$`Scénario07_50%`[[42]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_cov, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

# Biais median

dat <- data.frame(x_biais = c(comparaison04$`Scénario07_10%`[[1]], comparaison04$`Scénario07_50%`[[1]], comparaison04$`Scénario07_50%`[[1]],
                              comparaison04$`Scénario07_10%`[[7]], comparaison04$`Scénario07_50%`[[7]], comparaison04$`Scénario07_50%`[[7]],
                              comparaison04$`Scénario07_10%`[[13]], comparaison04$`Scénario07_50%`[[13]], comparaison04$`Scénario07_50%`[[13]],
                              comparaison04$`Scénario07_10%`[[19]], comparaison04$`Scénario07_50%`[[19]], comparaison04$`Scénario07_50%`[[19]],
                              comparaison04$`Scénario07_10%`[[25]], comparaison04$`Scénario07_50%`[[25]], comparaison04$`Scénario07_50%`[[25]],
                              comparaison04$`Scénario07_10%`[[31]], comparaison04$`Scénario07_50%`[[31]], comparaison04$`Scénario07_50%`[[31]],
                              comparaison04$`Scénario07_10%`[[37]], comparaison04$`Scénario07_50%`[[37]], comparaison04$`Scénario07_50%`[[37]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

################################################################################
################################################################################
# Scénarios avec différentes valeurs de prévalence de W2 et un taux de co-infection entre 30% et 40%
# nsim = 1000, Scénarios 08 prévalence de W2  ~ 5%, 10% et 15%

## Chargement des résultats bruts

## Taille d'échantillon n = 1000

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario08_5%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario08_10%.RData")

## Taille d'échantillon n = 5000

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario08_5%_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario08_10%_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario08_20%_5000.RData")

# Tableau de comparaison

comparaison05 <- data.frame(matrix(ncol = 8, 
                                   nrow = 42))

colnames(comparaison05) <- c(
  "Estimation", 
  "Performance",
  "Scénario08_5%",
  "Scénario08_10%",
  "Scénario08_20%",
  "Scénario08_5%_5000",
  "Scénario08_10%_5000",
  "Scénario08_20%_5000"
)

comparaison05$Estimation <- c("RegLog", "-", "-", "-", "-","-", "IPW", "-", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-", "-",
                              "TNDDR_GLM", "-", "-", "-", "-", "-")

comparaison05$Performance <- rep(c("Biais_med", "Biais_moy", "Variance", "MSE", "Précision", "%Cov"), 7)

comparaison05$`Scénario08_5%` <- c(
  
  (median(resultats_8_1$RRc) - mean(l_vraiRRc_8_1)), Tab01_8_1$Autres, Tab01_8_1$`Erreur de Monte Carlo`[4],
  (median(resultats2_8_1$RRm) - mean(l_vraiRRm_8_1)), Tab02_8_1$Autres, Tab02_8_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_1$RRm_RF, na.rm = TRUE) - mean(l_vraiRRm_8_1)), Tab03_8_1$Autres[1:4], Tab03_8_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_1$RRm_Lasso, na.rm = TRUE) - mean(l_vraiRRm_8_1)),Tab03_8_1$Autres[5:8], Tab03_8_1$`Erreur de Monte Carlo`[8],
  (median(resultats3_8_1$RRm_Mars) - mean(l_vraiRRm_8_1)), Tab03_8_1$Autres[9:12], Tab03_8_1$`Erreur de Monte Carlo`[12],
  (median(resultats3_8_1$RRm_RN) - mean(l_vraiRRm_8_1)), Tab03_8_1$Autres[13:16], Tab03_8_1$`Erreur de Monte Carlo`[16],
  (median(resultats3_8_1$RRm_GLM) - mean(l_vraiRRm_8_1)), Tab03_8_1$Autres[17:20], Tab03_8_1$`Erreur de Monte Carlo`[20]
  
)

comparaison05$`Scénario08_10%` <- c(
  
  (median(resultats_8_2$RRc) - mean(l_vraiRRc_8_2)), Tab01_8_2$Autres, Tab01_8_2$`Erreur de Monte Carlo`[4],
  (median(resultats2_8_2$RRm, na.rm = TRUE) - mean(l_vraiRRm_8_2)), Tab02_8_2$Autres, Tab02_8_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_2$RRm_RF, na.rm = TRUE) - mean(l_vraiRRm_8_2)), Tab03_8_2$Autres[1:4], Tab03_8_2$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_2$RRm_Lasso, na.rm = TRUE) - mean(l_vraiRRm_8_2)),Tab03_8_2$Autres[5:8], Tab03_8_2$`Erreur de Monte Carlo`[8],
  (median(resultats3_8_2$RRm_Mars) - mean(l_vraiRRm_8_2)), Tab03_8_2$Autres[9:12], Tab03_8_2$`Erreur de Monte Carlo`[12],
  (median(resultats3_8_2$RRm_RN) - mean(l_vraiRRm_8_2)), Tab03_8_2$Autres[13:16], Tab03_8_2$`Erreur de Monte Carlo`[16],
  (median(resultats3_8_2$RRm_GLM) - mean(l_vraiRRm_8_2)), Tab03_8_2$Autres[17:20], Tab03_8_2$`Erreur de Monte Carlo`[20]
  
)


comparaison05$`Scénario08_5%_5000` <- c(
  
  (median(resultats_8_1_5000$RRc) - mean(l_vraiRRc_8_1_5000)), Tab01_8_1_5000$Autres, Tab01_8_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_8_1_5000$RRm) - mean(l_vraiRRm_8_1_5000)), Tab02_8_1_5000$Autres, Tab02_8_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_1_5000$RRm_RF) - mean(l_vraiRRm_8_1_5000)), Tab03_8_1_5000$Autres[1:4], Tab03_8_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_1_5000$RRm_Lasso) - mean(l_vraiRRm_8_1_5000)),Tab03_8_1_5000$Autres[5:8], Tab03_8_1_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_8_1_5000$RRm_Mars) - mean(l_vraiRRm_8_1_5000)), Tab03_8_1_5000$Autres[9:12], Tab03_8_1_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_8_1_5000$RRm_RN) - mean(l_vraiRRm_8_1_5000)), Tab03_8_1_5000$Autres[13:16], Tab03_8_1_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_8_1_5000$RRm_GLM) - mean(l_vraiRRm_8_1_5000)), Tab03_8_1_5000$Autres[17:20], Tab03_8_1_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison05$`Scénario08_10%_5000` <- c(
  
  (median(resultats_8_2_5000$RRc) - mean(l_vraiRRc_8_2_5000)), Tab01_8_2_5000$Autres, Tab01_8_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_8_2_5000$RRm) - mean(l_vraiRRm_8_2_5000)), Tab02_8_2_5000$Autres, Tab02_8_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_2_5000$RRm_RF) - mean(l_vraiRRm_8_2_5000)), Tab03_8_2_5000$Autres[1:4], Tab03_8_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_2_5000$RRm_Lasso) - mean(l_vraiRRm_8_2_5000)),Tab03_8_2_5000$Autres[5:8], Tab03_8_2_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_8_2_5000$RRm_Mars) - mean(l_vraiRRm_8_2_5000)), Tab03_8_2_5000$Autres[9:12], Tab03_8_2_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_8_2_5000$RRm_RN) - mean(l_vraiRRm_8_2_5000)), Tab03_8_2_5000$Autres[13:16], Tab03_8_2_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_8_2_5000$RRm_GLM) - mean(l_vraiRRm_8_2_5000)), Tab03_8_2_5000$Autres[17:20], Tab03_8_2_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison05$`Scénario08_20%_5000` <- c(
  
  (median(resultats_8_3_5000$RRc) - mean(l_vraiRRc_8_3_5000)), Tab01_8_3_5000$Autres, Tab01_8_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_8_3_5000$RRm) - mean(l_vraiRRm_8_3_5000)), Tab02_8_3_5000$Autres, Tab02_8_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_3_5000$RRm_RF, na.rm = TRUE) - mean(l_vraiRRm_8_3_5000)), Tab03_8_3_5000$Autres[1:4], Tab03_8_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_8_3_5000$RRm_Lasso, na.rm = TRUE) - mean(l_vraiRRm_8_3_5000)),Tab03_8_3_5000$Autres[5:8], Tab03_8_3_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_8_3_5000$RRm_Mars) - mean(l_vraiRRm_8_3_5000)), Tab03_8_3_5000$Autres[9:12], Tab03_8_3_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_8_3_5000$RRm_RN) - mean(l_vraiRRm_8_3_5000)), Tab03_8_3_5000$Autres[13:16], Tab03_8_3_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_8_3_5000$RRm_GLM) - mean(l_vraiRRm_8_3_5000)), Tab03_8_3_5000$Autres[17:20], Tab03_8_3_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison05$`Scénario08_5%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison05$`Scénario08_5%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison05$`Scénario08_10%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison05$`Scénario08_10%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

comparaison05$`Scénario08_5%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison05$`Scénario08_5%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison05$`Scénario08_10%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison05$`Scénario08_10%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison05$`Scénario08_20%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison05$`Scénario08_20%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison05)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison05$`Scénario08_5%`[[2]], comparaison05$`Scénario08_10%`[[2]], comparaison05$`Scénario08_20%`[[2]],
                              comparaison05$`Scénario08_5%_5000`[[2]], comparaison05$`Scénario08_10%_5000`[[2]], comparaison05$`Scénario08_20%_5000`[[2]],
                              
                              comparaison05$`Scénario08_5%`[[8]], comparaison05$`Scénario08_10%`[[8]], comparaison05$`Scénario08_20%`[[8]],
                              comparaison05$`Scénario08_5%_5000`[[8]], comparaison05$`Scénario08_10%_5000`[[8]], comparaison05$`Scénario08_20%_5000`[[8]],
                              
                              comparaison05$`Scénario08_5%`[[14]], comparaison05$`Scénario08_10%`[[14]], comparaison05$`Scénario08_20%`[[14]],
                              comparaison05$`Scénario08_5%_5000`[[14]], comparaison05$`Scénario08_10%_5000`[[14]], comparaison05$`Scénario08_20%_5000`[[14]],
                              
                              comparaison05$`Scénario08_5%`[[20]], comparaison05$`Scénario08_10%`[[20]], comparaison05$`Scénario08_20%`[[20]],
                              comparaison05$`Scénario08_5%_5000`[[20]], comparaison05$`Scénario08_10%_5000`[[20]], comparaison05$`Scénario08_20%_5000`[[20]],
                              
                              comparaison05$`Scénario08_5%`[[26]], comparaison05$`Scénario08_10%`[[26]], comparaison05$`Scénario08_20%`[[26]],
                              comparaison05$`Scénario08_5%_5000`[[26]], comparaison05$`Scénario08_10%_5000`[[26]], comparaison05$`Scénario08_20%_5000`[[26]],
                              
                              comparaison05$`Scénario08_5%`[[32]], comparaison05$`Scénario08_10%`[[32]], comparaison05$`Scénario08_20%`[[32]],
                              comparaison05$`Scénario08_5%_5000`[[32]], comparaison05$`Scénario08_10%_5000`[[32]], comparaison05$`Scénario08_20%_5000`[[32]],
                              
                              comparaison05$`Scénario08_5%`[[38]], comparaison05$`Scénario08_10%`[[38]], comparaison05$`Scénario08_20%`[[38]],
                              comparaison05$`Scénario08_5%_5000`[[38]], comparaison05$`Scénario08_10%_5000`[[38]], comparaison05$`Scénario08_20%_5000`[[38]]),
                  
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 6),
                  Taille = rep(c("1000", "5000"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

# Couverture

dat <- data.frame(x_cov = c(comparaison03$`scénario08_5%`[[6]], comparaison03$`scénario08_10%`[[6]], comparaison03$`scénario08_15%`[[6]],
                            comparaison03$`scénario08_5%`[[12]], comparaison03$`scénario08_10%`[[12]], comparaison03$`scénario08_15%`[[12]],
                            comparaison03$`scénario08_5%`[[18]], comparaison03$`scénario08_10%`[[18]], comparaison03$`scénario08_15%`[[18]],
                            comparaison03$`scénario08_5%`[[24]], comparaison03$`scénario08_10%`[[24]], comparaison03$`scénario08_15%`[[24]],
                            comparaison03$`scénario08_5%`[[30]], comparaison03$`scénario08_10%`[[30]], comparaison03$`scénario08_15%`[[30]],
                            comparaison03$`scénario08_5%`[[36]], comparaison03$`scénario08_10%`[[36]], comparaison03$`scénario08_15%`[[36]],
                            comparaison03$`scénario08_5%`[[42]], comparaison03$`scénario08_10%`[[42]], comparaison03$`scénario08_15%`[[42]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_cov, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

# Biais median

dat <- data.frame(x_biais = c(comparaison03$`scénario08_5%`[[1]], comparaison03$`scénario08_10%`[[1]], comparaison03$`scénario08_15%`[[1]],
                              comparaison03$`scénario08_5%`[[7]], comparaison03$`scénario08_10%`[[7]], comparaison03$`scénario08_15%`[[7]],
                              comparaison03$`scénario08_5%`[[13]], comparaison03$`scénario08_10%`[[13]], comparaison03$`scénario08_15%`[[13]],
                              comparaison03$`scénario08_5%`[[19]], comparaison03$`scénario08_10%`[[19]], comparaison03$`scénario08_15%`[[19]],
                              comparaison03$`scénario08_5%`[[25]], comparaison03$`scénario08_10%`[[25]], comparaison03$`scénario08_15%`[[25]],
                              comparaison03$`scénario08_5%`[[31]], comparaison03$`scénario08_10%`[[31]], comparaison03$`scénario08_15%`[[31]],
                              comparaison03$`scénario08_5%`[[37]], comparaison03$`scénario08_10%`[[37]], comparaison03$`scénario08_15%`[[37]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

################################################################################
################################################################################
# Scénarios avec différentes valeurs de prévalence de W1 et un taux de co-infection entre 30% et 40%
# nsim = 1000, Scénarios 09 prévalence de W1  ~ 10%, 15% et 30%

## Chargement des résultats bruts

## Taille d'échantillon n = 1000

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario09_10%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario09_15%.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario09_20%.RData")

## Taille d'échantillon n = 5000

load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario09_10%_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario09_15%_5000.RData")
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'été 2026/Activité de recherche_04/Etudes-de-simulation/Résultats bruts/Résultats_Scénario09_20%_5000.RData")


# Tableau de comparaison

comparaison06 <- data.frame(matrix(ncol = 8, 
                                   nrow = 42))

colnames(comparaison06) <- c(
  "Estimation", 
  "Performance",
  "Scénario09_5%",
  "Scénario09_10%",
  "Scénario09_20%",
  "Scénario09_5%_5000",
  "Scénario09_10%_5000",
  "Scénario09_20%_5000"
)

comparaison06$Estimation <- c("RegLog", "-", "-", "-", "-","-", "IPW", "-", "-", "-", "-", "-", 
                              "TNDDR_RF", "-", "-", "-", "-", "-", "TNDDR_Lasso", "-", "-", "-", "-", "-",
                              "TNDDR_Mars", "-", "-", "-", "-", "-", "TNDDR_RN", "-", "-", "-", "-", "-",
                              "TNDDR_GLM", "-", "-", "-", "-", "-")

comparaison06$Performance <- rep(c("Biais_med", "Biais_moy", "Variance", "MSE", "Précision", "%Cov"), 7)

comparaison06$`Scénario09_5%` <- c(
  
  (median(resultats_9_1$RRc) - mean(l_vraiRRc_9_1)), Tab01_9_1$Autres, Tab01_9_1$`Erreur de Monte Carlo`[4],
  (median(resultats2_9_1$RRm) - mean(l_vraiRRm_9_1)), Tab02_9_1$Autres, Tab02_9_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_1$RRm_RF) - mean(l_vraiRRm_9_1)), Tab03_9_1$Autres[1:4], Tab03_9_1$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_1$RRm_Lasso) - mean(l_vraiRRm_9_1)),Tab03_9_1$Autres[5:8], Tab03_9_1$`Erreur de Monte Carlo`[8],
  (median(resultats3_9_1$RRm_Mars) - mean(l_vraiRRm_9_1)), Tab03_9_1$Autres[9:12], Tab03_9_1$`Erreur de Monte Carlo`[12],
  (median(resultats3_9_1$RRm_RN) - mean(l_vraiRRm_9_1)), Tab03_9_1$Autres[13:16], Tab03_9_1$`Erreur de Monte Carlo`[16],
  (median(resultats3_9_1$RRm_GLM) - mean(l_vraiRRm_9_1)), Tab03_9_1$Autres[17:20], Tab03_9_1$`Erreur de Monte Carlo`[20]
  
)

comparaison06$`Scénario09_10%` <- c(
  
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

comparaison06$`Scénario09_5%_5000` <- c(
  
  (median(resultats_9_1_5000$RRc) - mean(l_vraiRRc_9_1_5000)), Tab01_9_1_5000$Autres, Tab01_9_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_9_1_5000$RRm) - mean(l_vraiRRm_9_1_5000)), Tab02_9_1_5000$Autres, Tab02_9_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_1_5000$RRm_RF) - mean(l_vraiRRm_9_1_5000)), Tab03_9_1_5000$Autres[1:4], Tab03_9_1_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_1_5000$RRm_Lasso) - mean(l_vraiRRm_9_1_5000)),Tab03_9_1_5000$Autres[5:8], Tab03_9_1_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_9_1_5000$RRm_Mars) - mean(l_vraiRRm_9_1_5000)), Tab03_9_1_5000$Autres[9:12], Tab03_9_1_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_9_1_5000$RRm_RN) - mean(l_vraiRRm_9_1_5000)), Tab03_9_1_5000$Autres[13:16], Tab03_9_1_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_9_1_5000$RRm_GLM) - mean(l_vraiRRm_9_1_5000)), Tab03_9_1_5000$Autres[17:20], Tab03_9_1_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison06$`Scénario09_10%_5000` <- c(
  
  (median(resultats_9_2_5000$RRc) - mean(l_vraiRRc_9_2_5000)), Tab01_9_2_5000$Autres, Tab01_9_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_9_2_5000$RRm) - mean(l_vraiRRm_9_2_5000)), Tab02_9_2_5000$Autres, Tab02_9_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_2_5000$RRm_RF) - mean(l_vraiRRm_9_2_5000)), Tab03_9_2_5000$Autres[1:4], Tab03_9_2_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_2_5000$RRm_Lasso) - mean(l_vraiRRm_9_2_5000)),Tab03_9_2_5000$Autres[5:8], Tab03_9_2_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_9_2_5000$RRm_Mars) - mean(l_vraiRRm_9_2_5000)), Tab03_9_2_5000$Autres[9:12], Tab03_9_2_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_9_2_5000$RRm_RN) - mean(l_vraiRRm_9_2_5000)), Tab03_9_2_5000$Autres[13:16], Tab03_9_2_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_9_2_5000$RRm_GLM) - mean(l_vraiRRm_9_2_5000)), Tab03_9_2_5000$Autres[17:20], Tab03_9_2_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison06$`Scénario09_20%_5000` <- c(
  
  (median(resultats_9_3_5000$RRc) - mean(l_vraiRRc_9_3_5000)), Tab01_9_3_5000$Autres, Tab01_9_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats2_9_3_5000$RRm) - mean(l_vraiRRm_9_3_5000)), Tab02_9_3_5000$Autres, Tab02_9_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_3_5000$RRm_RF) - mean(l_vraiRRm_9_3_5000)), Tab03_9_3_5000$Autres[1:4], Tab03_9_3_5000$`Erreur de Monte Carlo`[4],
  (median(resultats3_9_3_5000$RRm_Lasso) - mean(l_vraiRRm_9_3_5000)),Tab03_9_3_5000$Autres[5:8], Tab03_9_3_5000$`Erreur de Monte Carlo`[8],
  (median(resultats3_9_3_5000$RRm_Mars) - mean(l_vraiRRm_9_3_5000)), Tab03_9_3_5000$Autres[9:12], Tab03_9_3_5000$`Erreur de Monte Carlo`[12],
  (median(resultats3_9_3_5000$RRm_RN) - mean(l_vraiRRm_9_3_5000)), Tab03_9_3_5000$Autres[13:16], Tab03_9_3_5000$`Erreur de Monte Carlo`[16],
  (median(resultats3_9_3_5000$RRm_GLM) - mean(l_vraiRRm_9_3_5000)), Tab03_9_3_5000$Autres[17:20], Tab03_9_3_5000$`Erreur de Monte Carlo`[20]
  
)

comparaison06$`Scénario09_5%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison06$`Scénario09_5%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison06$`Scénario09_10%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison06$`Scénario09_10%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison06$`Scénario09_20%`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison06$`Scénario09_20%`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

comparaison06$`Scénario09_5%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison06$`Scénario09_5%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison06$`Scénario09_10%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison06$`Scénario09_10%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])
comparaison06$`Scénario09_20%_5000`[-c(1, 7, 13, 19, 25, 31, 37)] <- sapply(comparaison06$`Scénario09_20%_5000`[-c(1, 7, 13, 19, 25, 31, 37)], FUN = function(x)x[[2]])

View(comparaison06)

####################### Représentation graphique ###############################

# Biais moyen

dat <- data.frame(x_biais = c(comparaison06$`Scénario09_5%`[[2]], comparaison06$`Scénario09_10%`[[2]], comparaison06$`Scénario09_20%`[[2]],
                              comparaison06$`Scénario09_5%_5000`[[2]], comparaison06$`Scénario09_10%_5000`[[2]], comparaison06$`Scénario09_20%_5000`[[2]],
                              
                              comparaison06$`Scénario09_5%`[[8]], comparaison06$`Scénario09_10%`[[8]], comparaison06$`Scénario09_20%`[[8]],
                              comparaison06$`Scénario09_5%_5000`[[8]], comparaison06$`Scénario09_10%_5000`[[8]], comparaison06$`Scénario09_20%_5000`[[8]],
                              
                              comparaison06$`Scénario09_5%`[[14]], comparaison06$`Scénario09_10%`[[14]], comparaison06$`Scénario09_20%`[[14]],
                              comparaison06$`Scénario09_5%_5000`[[14]], comparaison06$`Scénario09_10%_5000`[[14]], comparaison06$`Scénario09_20%_5000`[[14]],
                              
                              comparaison06$`Scénario09_5%`[[20]], comparaison06$`Scénario09_10%`[[20]], comparaison06$`Scénario09_20%`[[20]],
                              comparaison06$`Scénario09_5%_5000`[[20]], comparaison06$`Scénario09_10%_5000`[[20]], comparaison06$`Scénario09_20%_5000`[[20]],
                              
                              comparaison06$`Scénario09_5%`[[26]], comparaison06$`Scénario09_10%`[[26]], comparaison06$`Scénario09_20%`[[26]],
                              comparaison06$`Scénario09_5%_5000`[[26]], comparaison06$`Scénario09_10%_5000`[[26]], comparaison06$`Scénario09_20%_5000`[[26]],
                              
                              comparaison06$`Scénario09_5%`[[32]], comparaison06$`Scénario09_10%`[[32]], comparaison06$`Scénario09_20%`[[32]],
                              comparaison06$`Scénario09_5%_5000`[[32]], comparaison06$`Scénario09_10%_5000`[[32]], comparaison06$`Scénario09_20%_5000`[[32]],
                              
                              comparaison06$`Scénario09_5%`[[38]], comparaison06$`Scénario09_10%`[[38]], comparaison06$`Scénario09_20%`[[38]],
                              comparaison06$`Scénario09_5%_5000`[[38]], comparaison06$`Scénario09_10%_5000`[[38]], comparaison06$`Scénario09_20%_5000`[[38]]),
                  
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 6),
                  Taille = rep(c("1000", "5000"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

# Couverture

dat <- data.frame(x_cov = c(comparaison06$`Scénario09_5%`[[6]], comparaison06$`Scénario09_10%`[[6]], comparaison06$`Scénario09_15%`[[6]],
                            comparaison06$`Scénario09_5%`[[12]], comparaison06$`Scénario09_10%`[[12]], comparaison06$`Scénario09_15%`[[12]],
                            comparaison06$`Scénario09_5%`[[18]], comparaison06$`Scénario09_10%`[[18]], comparaison06$`Scénario09_15%`[[18]],
                            comparaison06$`Scénario09_5%`[[24]], comparaison06$`Scénario09_10%`[[24]], comparaison06$`Scénario09_15%`[[24]],
                            comparaison06$`Scénario09_5%`[[30]], comparaison06$`Scénario09_10%`[[30]], comparaison06$`Scénario09_15%`[[30]],
                            comparaison06$`Scénario09_5%`[[36]], comparaison06$`Scénario09_10%`[[36]], comparaison06$`Scénario09_15%`[[36]],
                            comparaison06$`Scénario09_5%`[[42]], comparaison06$`Scénario09_10%`[[42]], comparaison06$`Scénario09_15%`[[42]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_cov, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 

# Biais median

dat <- data.frame(x_biais = c(comparaison06$`Scénario09_5%`[[1]], comparaison06$`Scénario09_10%`[[1]], comparaison06$`Scénario09_15%`[[1]],
                              comparaison06$`Scénario09_5%`[[7]], comparaison06$`Scénario09_10%`[[7]], comparaison06$`Scénario09_15%`[[7]],
                              comparaison06$`Scénario09_5%`[[13]], comparaison06$`Scénario09_10%`[[13]], comparaison06$`Scénario09_15%`[[13]],
                              comparaison06$`Scénario09_5%`[[19]], comparaison06$`Scénario09_10%`[[19]], comparaison06$`Scénario09_15%`[[19]],
                              comparaison06$`Scénario09_5%`[[25]], comparaison06$`Scénario09_10%`[[25]], comparaison06$`Scénario09_15%`[[25]],
                              comparaison06$`Scénario09_5%`[[31]], comparaison06$`Scénario09_10%`[[31]], comparaison06$`Scénario09_15%`[[31]],
                              comparaison06$`Scénario09_5%`[[37]], comparaison06$`Scénario09_10%`[[37]], comparaison06$`Scénario09_15%`[[37]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN", "TNDDR_GLM"), each = 3),
                  Scénario = c(1, 2, 3))


dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


ggplot(dat, aes(x = Scénario, y = x_biais, color = Méthode,
                group = Méthode, shape = Méthode)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  
  labs(x = "Scénario", y = "Biais") +
  
  facet_wrap(~Taille)

theme_bw() 