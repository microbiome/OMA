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
tse <- transformSamples(tse, method = "relabundance", pseudocount = 1)
tse <- transformSamples(tse, assay_name = "relabundance", method = "clr")

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
##      Mixed    14  10
##      Veg       9  14
##                                         
##                Accuracy : 0.596         
##                  95% CI : (0.443, 0.736)
##     No Information Rate : 0.511         
##     P-Value [Acc > NIR] : 0.154         
##                                         
##                   Kappa : 0.192         
##                                         
##  Mcnemar's Test P-Value : 1.000         
##                                         
##             Sensitivity : 0.609         
##             Specificity : 0.583         
##          Pos Pred Value : 0.583         
##          Neg Pred Value : 0.609         
##              Prevalence : 0.489         
##          Detection Rate : 0.298         
##    Detection Prevalence : 0.511         
##       Balanced Accuracy : 0.596         
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
res <- evalm(model)
```

![](40_machine_learning_files/figure-latex/super5-1.pdf)<!-- --> ![](40_machine_learning_files/figure-latex/super5-2.pdf)<!-- --> ![](40_machine_learning_files/figure-latex/super5-3.pdf)<!-- --> ![](40_machine_learning_files/figure-latex/super5-4.pdf)<!-- --> 

```r
# Use patchwork to plot ROC and precision-recall curve side-by-side
library(patchwork)
res$roc + res$proc + 
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')
```

![](40_machine_learning_files/figure-latex/super5-5.pdf)<!-- --> 

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
 [5] ggplot2_3.3.6                  mikropml_1.3.0                
 [7] mia_1.5.16                     MultiAssayExperiment_1.22.0   
 [9] TreeSummarizedExperiment_2.1.4 Biostrings_2.64.1             
[11] XVector_0.36.0                 SingleCellExperiment_1.18.1   
[13] SummarizedExperiment_1.26.1    Biobase_2.56.0                
[15] GenomicRanges_1.48.0           GenomeInfoDb_1.32.4           
[17] IRanges_2.30.1                 S4Vectors_0.34.0              
[19] BiocGenerics_0.42.0            MatrixGenerics_1.8.1          
[21] matrixStats_0.62.0-9003        BiocStyle_2.24.0              
[23] rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] plyr_1.8.7                  lazyeval_0.2.2             
  [3] splines_4.2.1               BiocParallel_1.30.3        
  [5] listenv_0.8.0               scater_1.24.0              
  [7] digest_0.6.29               foreach_1.5.2              
  [9] yulab.utils_0.0.5           htmltools_0.5.3            
 [11] viridis_0.6.2               fansi_1.0.3                
 [13] magrittr_2.0.3              memoise_2.0.1              
 [15] MLmetrics_1.1.1             ScaledMatrix_1.4.1         
 [17] cluster_2.1.4               ROCR_1.0-11                
 [19] DECIPHER_2.24.0             recipes_1.0.1              
 [21] globals_0.16.1              gower_1.0.0                
 [23] hardhat_1.2.0               colorspace_2.0-3           
 [25] blob_1.2.3                  ggrepel_0.9.1              
 [27] xfun_0.33                   dplyr_1.0.10               
 [29] crayon_1.5.2                RCurl_1.98-1.9             
 [31] jsonlite_1.8.2              graph_1.74.0               
 [33] survival_3.4-0              iterators_1.0.14           
 [35] ape_5.6-2                   glue_1.6.2                 
 [37] gtable_0.3.1                ipred_0.9-13               
 [39] zlibbioc_1.42.0             DelayedArray_0.22.0        
 [41] kernlab_0.9-31              BiocSingular_1.12.0        
 [43] shape_1.4.6                 future.apply_1.9.1         
 [45] scales_1.2.1                DBI_1.1.3                  
 [47] Rcpp_1.0.9                  viridisLite_0.4.1          
 [49] decontam_1.16.0             tidytree_0.4.1             
 [51] proxy_0.4-27                bit_4.0.4                  
 [53] rsvd_1.0.5                  lava_1.6.10                
 [55] prodlim_2019.11.13          glmnet_4.1-4               
 [57] dir.expiry_1.4.0            ellipsis_0.3.2             
 [59] farver_2.1.1                pkgconfig_2.0.3            
 [61] XML_3.99-0.11               scuttle_1.6.3              
 [63] nnet_7.3-18                 CodeDepends_0.6.5          
 [65] utf8_1.2.2                  labeling_0.4.2             
 [67] tidyselect_1.1.2            rlang_1.0.6                
 [69] reshape2_1.4.4              munsell_0.5.0              
 [71] tools_4.2.1                 cachem_1.0.6               
 [73] xgboost_1.6.0.1             cli_3.4.1                  
 [75] DirichletMultinomial_1.38.0 generics_0.1.3             
 [77] RSQLite_2.2.18              evaluate_0.16              
 [79] stringr_1.4.1               fastmap_1.1.0              
 [81] yaml_2.3.5                  ModelMetrics_1.2.2.2       
 [83] knitr_1.40                  bit64_4.0.5                
 [85] randomForest_4.7-1.1        purrr_0.3.4                
 [87] future_1.28.0               nlme_3.1-159               
 [89] sparseMatrixStats_1.8.0     compiler_4.2.1             
 [91] beeswarm_0.4.0              filelock_1.0.2             
 [93] e1071_1.7-11                treeio_1.20.2              
 [95] tibble_3.1.8                stringi_1.7.8              
 [97] highr_0.9                   Matrix_1.5-1               
 [99] vegan_2.6-2                 permute_0.9-7              
[101] vctrs_0.4.2                 pillar_1.8.1               
[103] lifecycle_1.0.2             BiocManager_1.30.18        
[105] BiocNeighbors_1.14.0        data.table_1.14.2          
[107] bitops_1.0-7                irlba_2.3.5.1              
[109] R6_2.5.1                    bookdown_0.29              
[111] gridExtra_2.3               vipor_0.4.5                
[113] parallelly_1.32.1           codetools_0.2-18           
[115] MASS_7.3-58.1               assertthat_0.2.1           
[117] withr_2.5.0                 GenomeInfoDbData_1.2.8     
[119] mgcv_1.8-40                 parallel_4.2.1             
[121] grid_4.2.1                  rpart_4.1.16               
[123] beachmat_2.12.0             timeDate_4021.106          
[125] tidyr_1.2.1                 class_7.3-20               
[127] rmarkdown_2.16              DelayedMatrixStats_1.18.1  
[129] pROC_1.18.0                 lubridate_1.8.0            
[131] ggbeeswarm_0.6.0           
```
</div>

