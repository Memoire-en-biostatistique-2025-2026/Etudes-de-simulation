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
load("C:/Users/lenovo/Desktop/Université Laval 2025-2026/Session d'hiver 2026/Activité de recherche_03/Etudes-de-simulation/Résultats_Scénario03_40%.RData")

# Tableau de comparaison 

comparaison01 <- data.frame(matrix(ncol = 6, 
                                   nrow = 30))

colnames(comparaison01) <- c(
  "Estimation", 
  "Performance",
  "Scénario01",
  "Scénario02",
  "Scénario03_20%",
  "Scénario03_40%"
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

comparaison01$`Scénario03_20%` <- c(
  
  Tab1_3$Autres, Tab1_3$`Erreur de Monte Carlo`[4],
  Tab2_3$Autres, Tab2_3$`Erreur de Monte Carlo`[4],
  Tab3_3$Autres[1:4], Tab3_3$`Erreur de Monte Carlo`[4],
  Tab3_3$Autres[5:8], Tab3_3$`Erreur de Monte Carlo`[8],
  Tab3_3$Autres[9:12], Tab3_3$`Erreur de Monte Carlo`[12],
  Tab3_3$Autres[13:16], Tab3_3$`Erreur de Monte Carlo`[16]
  
)

comparaison01$`Scénario03_40%` <- c(
  
  tab1_4$Autres, tab1_4$`Erreur de Monte Carlo`[4],
  tab2_4$Autres, tab2_4$`Erreur de Monte Carlo`[4],
  tab3_4$Autres[1:4], tab3_4$`Erreur de Monte Carlo`[4],
  tab3_4$Autres[5:8], tab3_4$`Erreur de Monte Carlo`[8],
  tab3_4$Autres[9:12], tab3_4$`Erreur de Monte Carlo`[12],
  tab3_4$Autres[13:16], tab3_4$`Erreur de Monte Carlo`[16]
  
)

# nsim = 1000, couverture vaccinale  ~ 33% : Scénarios 01, 02 et 03

comparaison01$Scénario01 <- sapply(comparaison01$Scénario01, FUN = function(x)x[[2]])
comparaison01$Scénario02 <- sapply(comparaison01$Scénario02, FUN = function(x)x[[2]][[1]])
comparaison01$`Scénario03_20%` <- sapply(comparaison01$`Scénario03_20%`, FUN = function(x)x[[2]][[1]])
comparaison01$`Scénario03_40%` <- sapply(comparaison01$`Scénario03_40%`, FUN = function(x)x[[2]][[1]])

View(comparaison01)

# nsim = 1000, coouverture vaccinale  ~ 33%, 50% et 70% : Scénarios 02 et 03

comparaison02 <- data.frame(matrix(ncol = 11, 
                                   nrow = 30))

colnames(comparaison02) <- c(
  "Estimation", 
  "Performance",
  "Scénario02_33%",
  "Scénario02_50%",
  "Scénario02_70%",
  "Scénario03_33%",
  "Scénario03_50%",
  "Scénario03_70%",
  "Scénario04_33%",
  "Scénario04_50%",
  "Scénario04_70%"
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

comparaison02$`Scénario04_33%` <- c(
  
  sapply(tab1_4$Autres, function(x) x[[2]][[1]]), tab1_4$`Erreur de Monte Carlo`[[4]][[2]][[1]],
  sapply(tab1_4$Autres, function(x) x[[2]][[1]]), tab2_4$`Erreur de Monte Carlo`[[4]][[2]][[1]],
  sapply(tab3_4$Autres[1:4],function(x) x[[2]][[1]]) , tab3_4$`Erreur de Monte Carlo`[[4]][[2]][[1]],
  sapply(tab3_4$Autres[5:8],function(x) x[[2]][[1]]), tab3_4$`Erreur de Monte Carlo`[[8]][[2]][[1]],
  sapply(tab3_4$Autres[9:12],function(x) x[[2]][[1]]), tab3_4$`Erreur de Monte Carlo`[[12]][[2]][[1]],
  sapply(tab3_4$Autres[13:16],function(x) x[[2]][[1]]), tab3_4$`Erreur de Monte Carlo`[[16]][[2]][[1]]
  
)

comparaison02$`Scénario04_50%` <- c(
  
  sapply(tab1_4$Autres, function(x) x[[2]][[2]]), tab1_4$`Erreur de Monte Carlo`[[4]][[2]][[2]],
  sapply(tab1_4$Autres, function(x) x[[2]][[2]]), tab2_4$`Erreur de Monte Carlo`[[4]][[2]][[2]],
  sapply(tab3_4$Autres[1:4],function(x) x[[2]][[2]]) , tab3_4$`Erreur de Monte Carlo`[[4]][[2]][[2]],
  sapply(tab3_4$Autres[5:8],function(x) x[[2]][[2]]), tab3_4$`Erreur de Monte Carlo`[[8]][[2]][[2]],
  sapply(tab3_4$Autres[9:12],function(x) x[[2]][[2]]), tab3_4$`Erreur de Monte Carlo`[[12]][[2]][[2]],
  sapply(tab3_4$Autres[13:16],function(x) x[[2]][[2]]), tab3_4$`Erreur de Monte Carlo`[[16]][[2]][[2]]
  
)

comparaison02$`Scénario04_70%` <- c(
  
  sapply(tab1_4$Autres, function(x) x[[2]][[3]]), tab1_4$`Erreur de Monte Carlo`[[4]][[2]][[3]],
  sapply(tab1_4$Autres, function(x) x[[2]][[3]]), tab2_4$`Erreur de Monte Carlo`[[4]][[2]][[3]],
  sapply(tab3_4$Autres[1:4],function(x) x[[2]][[3]]) , tab3_4$`Erreur de Monte Carlo`[[4]][[2]][[3]],
  sapply(tab3_4$Autres[5:8],function(x) x[[2]][[3]]), tab3_4$`Erreur de Monte Carlo`[[8]][[2]][[3]],
  sapply(tab3_4$Autres[9:12],function(x) x[[2]][[3]]), tab3_4$`Erreur de Monte Carlo`[[12]][[2]][[3]],
  sapply(tab3_4$Autres[13:16],function(x) x[[2]][[3]]), tab3_4$`Erreur de Monte Carlo`[[16]][[2]][[3]]
  
)

View(comparaison02)
################################################################################

dat <- data.frame(x_biais = c(comparaison01$Scénario01[[1]], comparaison01$Scénario02[[1]], comparaison01$Scénario03[[1]],
                   comparaison01$Scénario01[[6]], comparaison01$Scénario02[[6]], comparaison01$Scénario03[[6]],
                   comparaison01$Scénario01[[11]], comparaison01$Scénario02[[11]], comparaison01$Scénario03[[11]],
                  comparaison01$Scénario01[[16]], comparaison01$Scénario02[[16]], comparaison01$Scénario03[[16]],
                  comparaison01$Scénario01[[21]], comparaison01$Scénario02[[21]], comparaison01$Scénario03[[21]],
                  comparaison01$Scénario01[[26]], comparaison01$Scénario02[[26]], comparaison01$Scénario03[[26]]),
                  Méthode = rep(c("RegLog", "IPW", "TNDDR_RF", "TNDDR_Lasso", "TNDDR_Mars", "TNDDR_RN"), each = 3),
                  Scénario = c(1, 2, 3))
                               

dat$Scénario <- as.factor(dat$Scénario)
dat$Méthode<- as.factor(dat$Méthode)


x <- ggplot(dat,aes(x = Scénario, y = x_biais, fill = Méthode))
y <- x + geom_bar(position = position_dodge(), stat = "identity")

y + labs(y = "Biais")

x <- ggplot(dat,aes(x = Scénario, y = x_biais, fill = Méthode))
y <- x + geom_point()
y + labs(y = "Moyenne de la masse corporelle (g)")
 
ggplot(dat) +
  aes(x = Scénario, y = x_biais) +
  geom_col(fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(Méthode))
