## ----setup, echo=FALSE, results="asis"-----------------
library(rebook)
chapterPreamble()


## ----install-pkg, include = FALSE----------------------
if(!require(mikropml)){
  install.packages("mikropml")
}

if(!require(MLeval)){
  install.packages("MLeval")
}


## ----superML1------------------------------------------
library(mia)

# Load experimental data
data(peerj13075, package="mia")
tse <- peerj13075


## ----super2--------------------------------------------
# Agglomerate data
tse <- mergeFeaturesByRank(tse, rank = "order")

# Apply CLR transform
tse <- transformAssay(tse, assay.type = "counts", method = "clr",
                       MARGIN="samples", pseudocount=1)

# Get assay
assay <- assay(tse, "clr")
# Transpose assay
assay <- t(assay)

# Convert into data.frame
df <- as.data.frame(assay)

# Add labels to assay
labels <- colData(tse)$Diet
labels <- as.factor(labels)
df$diet <- labels 

df[5, 5]


## ----super3--------------------------------------------
library(mikropml)

# Run random forest 
results <- run_ml(df, "rf", outcome_colname = "diet", 
                  kfold = 2, cv_times = 5, training_frac = 0.8)

# Print result
confusionMatrix(data = results$trained_model$finalModel$predicted, 
                reference = results$trained_model$finalModel$y)


## ----super4--------------------------------------------
# Set seed for reproducibility
set.seed(6358)

# Specify train control
train_control <- trainControl(method = "cv", number = 5,
                              classProbs = TRUE, 
                              savePredictions = "final",
                              allowParallel = TRUE)

# Specify hyperparameter tuning grid
tune_grid <- expand.grid(nrounds = c(50, 100, 200),
                         max_depth = c(6, 8, 10),
                         colsample_bytree = c(0.6, 0.8, 1),
                         eta = c(0.1, 0.3),
                         gamma = 0,
                         min_child_weight = c(3, 4, 5),
                         subsample = c(0.6, 0.8)
                         )

# Train the model, use LOOCV to evaluate performance
model <- train(x = assay, 
               y = labels, 
               method = "xgbTree",
               objective = "binary:logistic",
               trControl = train_control,
               tuneGrid = tune_grid,
               metric = "AUC",
               verbosity = 0
)



## ----super5--------------------------------------------
library(MLeval)

# Calculate different evaluation metrics
res <- evalm(model, showplots = FALSE)

# Use patchwork to plot ROC and precision-recall curve side-by-side
library(patchwork)
res$roc + res$proc + 
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')

