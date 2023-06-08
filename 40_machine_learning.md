# Machine learning {#machine_learning}

<script>
document.addEventListener("click", function (event) {
    if (event.target.classList.contains("rebook-collapse")) {
        event.target.classList.toggle("active");
        var content = event.target.nextElementSibling;
        if (content.style.display === "block") {
            content.style.display = "none";
        } else {
            content.style.display = "block";
        }
    }
})
</script>

<style>
.rebook-collapse {
  background-color: #eee;
  color: #444;
  cursor: pointer;
  padding: 18px;
  width: 100%;
  border: none;
  text-align: left;
  outline: none;
  font-size: 15px;
}

.rebook-content {
  padding: 0 18px;
  display: none;
  overflow: hidden;
  background-color: #f1f1f1;
}
</style>

Machine learning (ML) is a part of artificial intelligence. There are multiple
definitions, but "machine" refers to computation and "learning" to improving 
performance based on the data by finding patterns from it. Machine learning
includes wide variety of methods from simple statistical methods to more 
complex methods such as neural-networks. 

Machine learning can be divided into supervised and unsupervised machine learning.
Supervised ML is used to predict outcome based on the data. Unsupervised ML is used, 
for example, to reduce dimensionality (e.g. PCA) and to find clusters from the 
data (e.g., k-means clustering). 


## Supervised machine learning



"Supervised" means that the training data is introduced before. The training data
contains labels (e.g., patient status), and the model is fitted based on the 
training data. After fitting, the model is utilized to predict labels of data whose 
labels are not known. 


```r
library(mia)

# Load experimental data
data(peerj13075, package="mia")
tse <- peerj13075
```

Let's first preprocess the data.


```r
# Agglomerate data
tse <- agglomerateByRank(tse, rank = "order")

# Apply CLR transform
tse <- transformCounts(tse, assay.type = "counts", method = "clr",
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
```

```
## [1] -0.4612
```

In the example below, we use [mikropml](https://journals.asm.org/doi/10.1128/mBio.00434-20)
package. We try to predict the diet type based on the data.


```r
library(mikropml)

# Run random forest 
results <- run_ml(df, "rf", outcome_colname = "diet", 
                  kfold = 2, cv_times = 5, training_frac = 0.8)

# Print result
confusionMatrix(data = results$trained_model$finalModel$predicted, 
                reference = results$trained_model$finalModel$y)
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction Mixed Veg
##      Mixed    16   6
##      Veg       7  18
##                                         
##                Accuracy : 0.723         
##                  95% CI : (0.574, 0.844)
##     No Information Rate : 0.511         
##     P-Value [Acc > NIR] : 0.00244       
##                                         
##                   Kappa : 0.446         
##                                         
##  Mcnemar's Test P-Value : 1.00000       
##                                         
##             Sensitivity : 0.696         
##             Specificity : 0.750         
##          Pos Pred Value : 0.727         
##          Neg Pred Value : 0.720         
##              Prevalence : 0.489         
##          Detection Rate : 0.340         
##    Detection Prevalence : 0.468         
##       Balanced Accuracy : 0.723         
##                                         
##        'Positive' Class : Mixed         
## 
```

mikropml offers easier interface to [caret](https://cran.r-project.org/web/packages/caret/index.html) 
package. However, we can also use it directly.

Let's use xgboost model which is another commonly used algorithm in bioinformatics.


```r
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
```

Let's create ROC curve which is a commonly used method in binary classification.
For unbalanced data, you might want to plot precision-recall curve. 


```r
library(MLeval)

# Calculate different evaluation metrics
res <- evalm(model, showplots = FALSE)

# Use patchwork to plot ROC and precision-recall curve side-by-side
library(patchwork)
res$roc + res$proc + 
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')
```

![](40_machine_learning_files/figure-latex/super5-1.pdf)<!-- --> 

## Unsupervised machine learning

"Unsupervised" means that the labels (e.g., patient status is not known), 
and patterns are learned based only the abundance table, for instance. 
Unsupervised ML is also known as a data mining where patterns are extracted 
from big datasets. 

For unsupervised machine learning, please refer to chapters that are listed below:

- Chapter \@ref(clustering)
- Chapter \@ref(community-similarity) 

## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.3.0 (2023-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] patchwork_1.1.2                caret_6.0-94                  
 [3] lattice_0.21-8                 ggplot2_3.4.2                 
 [5] mia_1.9.2                      MultiAssayExperiment_1.26.0   
 [7] TreeSummarizedExperiment_2.1.4 Biostrings_2.68.1             
 [9] XVector_0.40.0                 SingleCellExperiment_1.22.0   
[11] SummarizedExperiment_1.30.2    Biobase_2.60.0                
[13] GenomicRanges_1.52.0           GenomeInfoDb_1.36.0           
[15] IRanges_2.34.0                 S4Vectors_0.38.1              
[17] BiocGenerics_0.46.0            MatrixGenerics_1.12.0         
[19] matrixStats_1.0.0-9000         MLeval_0.3                    
[21] mikropml_1.6.0                 BiocStyle_2.28.0              
[23] rebook_1.10.1                 

loaded via a namespace (and not attached):
  [1] splines_4.3.0               bitops_1.0-7               
  [3] filelock_1.0.2              tibble_3.2.1               
  [5] CodeDepends_0.6.5           hardhat_1.3.0              
  [7] graph_1.78.0                pROC_1.18.2                
  [9] XML_3.99-0.14               rpart_4.1.19               
 [11] DirichletMultinomial_1.42.0 lifecycle_1.0.3            
 [13] globals_0.16.2              MASS_7.3-60                
 [15] magrittr_2.0.3              rmarkdown_2.22             
 [17] yaml_2.3.7                  DBI_1.1.3                  
 [19] lubridate_1.9.2             zlibbioc_1.46.0            
 [21] purrr_1.0.1                 RCurl_1.98-1.12            
 [23] yulab.utils_0.0.6           nnet_7.3-19                
 [25] ipred_0.9-14                lava_1.7.2.1               
 [27] GenomeInfoDbData_1.2.10     ggrepel_0.9.3              
 [29] irlba_2.3.5.1               listenv_0.9.0              
 [31] tidytree_0.4.2              vegan_2.6-4                
 [33] parallelly_1.36.0           permute_0.9-7              
 [35] DelayedMatrixStats_1.22.0   codetools_0.2-19           
 [37] DelayedArray_0.26.3         scuttle_1.10.1             
 [39] tidyselect_1.2.0            shape_1.4.6                
 [41] farver_2.1.1                ScaledMatrix_1.8.1         
 [43] viridis_0.6.3               jsonlite_1.8.5             
 [45] BiocNeighbors_1.18.0        e1071_1.7-13               
 [47] decontam_1.20.0             survival_3.5-5             
 [49] scater_1.28.0               iterators_1.0.14           
 [51] foreach_1.5.2               tools_4.3.0                
 [53] treeio_1.24.1               Rcpp_1.0.10                
 [55] glue_1.6.2                  prodlim_2023.03.31         
 [57] gridExtra_2.3               xfun_0.39                  
 [59] mgcv_1.8-42                 dplyr_1.1.2                
 [61] withr_2.5.0                 BiocManager_1.30.20        
 [63] fastmap_1.1.1               fansi_1.0.4                
 [65] digest_0.6.31               rsvd_1.0.5                 
 [67] timechange_0.2.0            R6_2.5.1                   
 [69] colorspace_2.1-0            RSQLite_2.3.1              
 [71] utf8_1.2.3                  tidyr_1.3.0                
 [73] generics_0.1.3              data.table_1.14.8          
 [75] recipes_1.0.6               DECIPHER_2.28.0            
 [77] class_7.3-22                S4Arrays_1.0.4             
 [79] ModelMetrics_1.2.2.2        pkgconfig_2.0.3            
 [81] gtable_0.3.3                timeDate_4022.108          
 [83] blob_1.2.4                  MLmetrics_1.1.1            
 [85] htmltools_0.5.5             bookdown_0.34              
 [87] scales_1.2.1                gower_1.0.1                
 [89] knitr_1.43                  rstudioapi_0.14            
 [91] reshape2_1.4.4              nlme_3.1-162               
 [93] proxy_0.4-27                cachem_1.0.8               
 [95] stringr_1.5.0               parallel_4.3.0             
 [97] vipor_0.4.5                 pillar_1.9.0               
 [99] grid_4.3.0                  vctrs_0.6.2                
[101] randomForest_4.7-1.1        BiocSingular_1.16.0        
[103] beachmat_2.16.0             cluster_2.1.4              
[105] beeswarm_0.4.0              evaluate_0.21              
[107] cli_3.6.1                   compiler_4.3.0             
[109] rlang_1.1.1                 crayon_1.5.2               
[111] future.apply_1.11.0         labeling_0.4.2             
[113] plyr_1.8.8                  ggbeeswarm_0.7.2           
[115] stringi_1.7.12              viridisLite_0.4.2          
[117] BiocParallel_1.34.2         munsell_0.5.0              
[119] lazyeval_0.2.2              glmnet_4.1-7               
[121] Matrix_1.5-4.1              dir.expiry_1.8.0           
[123] sparseMatrixStats_1.12.0    bit64_4.0.5                
[125] future_1.32.0               highr_0.10                 
[127] kernlab_0.9-32              ROCR_1.0-11                
[129] memoise_2.0.1               bit_4.0.5                  
[131] xgboost_1.7.5.1             ape_5.7-1                  
```
</div>

