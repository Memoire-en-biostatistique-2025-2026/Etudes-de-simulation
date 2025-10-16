# 1. Calculer les vraies valeurs (fait !)
# 2. Parametres de la simulation
#    - Nombre de replications (1000 ?)
#    - Liste de 1000 germes
# 3. Initialiser des objets pour contenir les resultats
#    - Coefficient de la regression logistique
#    - Erreur-type du coefficient
# 4. Faire une boucle de 1 a nombre de replications
# 5. Generation des donnees TND (n = ?)
# 6. Analyse avec regression logistique
# 7. Enregistrer les resultats
# 8. Analyser les resultats
#    - statistiques descriptives
#    - biais, variance, moyenne de l'erreur-type, couverture des IC
################################################################################
################################################################################

# 3. Initialiser des objets pour contenir les resultats
#    - Coefficient de la regression logistique
#    - Erreur-type du coefficient

coe_reg <- rep(NA, 1000)

err_reg <- rep(NA, 1000)

# 4. Faire une boucle de 1 a nombre de replications
# 5. Generation des donnees TND (n = 1000)

co_inf_para <- # Proportion des co-infections
  
  set.seed(1) # Pour avoir toujours les mêmes germes

seeds_list <- sample(1:1000000, size = 1000)

l_vraiRRc <- rep(NA, 1000)

for (i in seeds_list) {
  
  for (j in 1:1000) {
    
    TNDdat <- datagen(seed = i, ssize = 1000, co_inf_para = co_inf_para)
    
  }
}