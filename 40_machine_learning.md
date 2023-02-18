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
data(peerj13075)
(tse <- peerj13075)
```

```
## class: TreeSummarizedExperiment 
## dim: 674 58 
## metadata(0):
## assays(1): counts
## rownames(674): OTU1 OTU2 ... OTU2567 OTU2569
## rowData names(6): kingdom phylum ... family genus
## colnames(58): ID1 ID2 ... ID57 ID58
## colData names(5): Sample Geographical_location Gender Age Diet
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: NULL
## rowTree: NULL
## colLinks: NULL
## colTree: NULL
```

Let's first preprocess the data.


```r
# Agglomerate data
tse <- agglomerateByRank(tse, rank = "order")

# Apply CLR transform
assay(tse, "pseudo") <- assay(tse, "counts") + 1
tse <- transformCounts(tse, assay_name = "pseudo", method = "relabundance")
tse <- transformCounts(tse, assay_name = "relabundance", method = "clr")

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
if( !require("mikropml") ){
    install.packages("mikropml")
    library(mikropml)
}

# Run random forest 
results <- run_ml(df, "rf", outcome_colname = 'diet', 
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
##      Mixed    15   9
##      Veg       8  15
##                                         
##                Accuracy : 0.638         
##                  95% CI : (0.485, 0.773)
##     No Information Rate : 0.511         
##     P-Value [Acc > NIR] : 0.0536        
##                                         
##                   Kappa : 0.277         
##                                         
##  Mcnemar's Test P-Value : 1.0000        
##                                         
##             Sensitivity : 0.652         
##             Specificity : 0.625         
##          Pos Pred Value : 0.625         
##          Neg Pred Value : 0.652         
##              Prevalence : 0.489         
##          Detection Rate : 0.319         
##    Detection Prevalence : 0.511         
##       Balanced Accuracy : 0.639         
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
if( !require(MLeval) ){
    install.packages("MLeval")
    library(MLeval)
}
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
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] patchwork_1.1.2                MLeval_0.3                    
 [3] caret_6.0-93                   lattice_0.20-45               
 [5] ggplot2_3.4.1                  mikropml_1.5.0                
 [7] mia_1.7.5                      MultiAssayExperiment_1.24.0   
 [9] TreeSummarizedExperiment_2.1.4 Biostrings_2.66.0             
[11] XVector_0.38.0                 SingleCellExperiment_1.20.0   
[13] SummarizedExperiment_1.28.0    Biobase_2.58.0                
[15] GenomicRanges_1.50.2           GenomeInfoDb_1.34.9           
[17] IRanges_2.32.0                 S4Vectors_0.36.1              
[19] BiocGenerics_0.44.0            MatrixGenerics_1.10.0         
[21] matrixStats_0.63.0-9003        BiocStyle_2.24.0              
[23] rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] plyr_1.8.8                  lazyeval_0.2.2             
  [3] splines_4.2.1               BiocParallel_1.32.5        
  [5] listenv_0.9.0               scater_1.26.1              
  [7] digest_0.6.31               foreach_1.5.2              
  [9] yulab.utils_0.0.6           htmltools_0.5.4            
 [11] viridis_0.6.2               fansi_1.0.4                
 [13] magrittr_2.0.3              memoise_2.0.1              
 [15] MLmetrics_1.1.1             ScaledMatrix_1.6.0         
 [17] cluster_2.1.4               ROCR_1.0-11                
 [19] DECIPHER_2.26.0             recipes_1.0.4              
 [21] globals_0.16.2              gower_1.0.1                
 [23] hardhat_1.2.0               timechange_0.2.0           
 [25] colorspace_2.1-0            blob_1.2.3                 
 [27] ggrepel_0.9.3               xfun_0.37                  
 [29] dplyr_1.1.0                 crayon_1.5.2               
 [31] RCurl_1.98-1.10             jsonlite_1.8.4             
 [33] graph_1.74.0                survival_3.5-3             
 [35] iterators_1.0.14            ape_5.7                    
 [37] glue_1.6.2                  gtable_0.3.1               
 [39] ipred_0.9-13                zlibbioc_1.44.0            
 [41] DelayedArray_0.24.0         kernlab_0.9-32             
 [43] BiocSingular_1.14.0         shape_1.4.6                
 [45] future.apply_1.10.0         scales_1.2.1               
 [47] DBI_1.1.3                   Rcpp_1.0.10                
 [49] viridisLite_0.4.1           decontam_1.18.0            
 [51] tidytree_0.4.2              proxy_0.4-27               
 [53] bit_4.0.5                   rsvd_1.0.5                 
 [55] lava_1.7.1                  prodlim_2019.11.13         
 [57] glmnet_4.1-6                dir.expiry_1.4.0           
 [59] farver_2.1.1                pkgconfig_2.0.3            
 [61] XML_3.99-0.13               scuttle_1.8.4              
 [63] nnet_7.3-18                 CodeDepends_0.6.5          
 [65] utf8_1.2.3                  labeling_0.4.2             
 [67] tidyselect_1.2.0            rlang_1.0.6                
 [69] reshape2_1.4.4              munsell_0.5.0              
 [71] tools_4.2.1                 cachem_1.0.6               
 [73] xgboost_1.7.3.1             cli_3.6.0                  
 [75] DirichletMultinomial_1.40.0 generics_0.1.3             
 [77] RSQLite_2.2.20              evaluate_0.20              
 [79] stringr_1.5.0               fastmap_1.1.0              
 [81] yaml_2.3.7                  ModelMetrics_1.2.2.2       
 [83] knitr_1.42                  bit64_4.0.5                
 [85] randomForest_4.7-1.1        purrr_1.0.1                
 [87] future_1.31.0               nlme_3.1-162               
 [89] sparseMatrixStats_1.10.0    compiler_4.2.1             
 [91] beeswarm_0.4.0              filelock_1.0.2             
 [93] e1071_1.7-13                treeio_1.22.0              
 [95] tibble_3.1.8                stringi_1.7.12             
 [97] highr_0.10                  Matrix_1.5-3               
 [99] vegan_2.6-4                 permute_0.9-7              
[101] vctrs_0.5.2                 pillar_1.8.1               
[103] lifecycle_1.0.3             BiocManager_1.30.19        
[105] BiocNeighbors_1.16.0        data.table_1.14.6          
[107] bitops_1.0-7                irlba_2.3.5.1              
[109] R6_2.5.1                    bookdown_0.32              
[111] gridExtra_2.3               vipor_0.4.5                
[113] parallelly_1.34.0           codetools_0.2-19           
[115] MASS_7.3-58.2               withr_2.5.0                
[117] GenomeInfoDbData_1.2.9      mgcv_1.8-41                
[119] parallel_4.2.1              grid_4.2.1                 
[121] rpart_4.1.19                beachmat_2.14.0            
[123] timeDate_4022.108           tidyr_1.3.0                
[125] class_7.3-21                rmarkdown_2.20             
[127] DelayedMatrixStats_1.20.0   pROC_1.18.0                
[129] lubridate_1.9.2             ggbeeswarm_0.7.1           
```
</div>

