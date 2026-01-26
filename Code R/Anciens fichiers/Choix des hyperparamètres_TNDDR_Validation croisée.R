# Chargement des librairies nécessaires

library(caret)
library(mlbench)
library(ranger)


param_grid <- expand.grid(
                       mtry = 1,
                       splitrule = "extratrees",
                       min.node.size = 10:60
                       )

cv_scheme <- trainControl(method = "cv",
                          number = 5,
                          verboseIter = FALSE)
models <- list()

for (ntree in c(500, 800, 1000)){
  
  rf_model <- train(
                    x = TNDdat[, "C", drop = FALSE],
                    y = TNDdat[, 2],
                    method = "ranger",
                    trControl = cv_scheme,
                    tuneGrid = param_grid,
                    num.trees=ntree
                    )
  
  name <- paste0(ntree,"_tr_model")
  models[[name]] <- rf_model
}

models <- list()

for (ntree in c(50, 500, 800, 1000)){
  
  rf_model <- train(
    x = TNDdat_train1[, c("C", "V"), drop = FALSE],
    y = TNDdat_train1[, 3],
    method = "ranger",
    trControl = cv_scheme,
    tuneGrid = param_grid,
    num.trees=ntree
  )
  
  name <- paste0(ntree,"_tr_model")
  models[[name]] <- rf_model
}

