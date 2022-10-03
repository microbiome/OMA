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

head(df, 5)
```

```
##     Lactobacillales Bacillales Clostridiales Selenomonadales Acholeplasmatales
## ID1         3.71440      9.359         2.616          1.2295           -0.5623
## ID2         1.87823     10.238         2.616          0.8666           -0.5197
## ID3         0.07201      4.052         5.631         -0.6211           -0.6211
## ID4         1.25738      6.462         5.444          0.3411           -0.3521
## ID5         6.62219      8.639         1.841          0.2320           -0.4612
##     Burkholderiales Acidaminococcales Rhodospirillales Acidimicrobiales 
## ID1         -0.5623           -0.5623          -0.5623           -0.5623
## ID2         -0.5197            0.1735          -0.5197            0.1735
## ID3         -0.6211           -0.6211          -0.6211           -0.6211
## ID4         -0.3521           -0.3521          -0.3521           -0.3521
## ID5         -0.4612           -0.4612          -0.4612            0.9251
##     AcidithioBacillales Pseudomonadales      NA Actinomycetales
## ID1             -0.5623           5.986 -0.5623          0.1309
## ID2             -0.5197           3.144 -0.5197          0.1735
## ID3             -0.6211           5.454 -0.6211         -0.6211
## ID4             -0.3521           2.538 -0.3521         -0.3521
## ID5             -0.4612           4.765 -0.4612          0.2320
##     Micromonosporales Streptosporangiales Pseudonocardiales Actinopolysporales
## ID1           -0.5623            -0.56226            0.1309            -0.5623
## ID2            0.1735            -0.51966            0.1735            -0.5197
## ID3           -0.6211             0.07201            0.7652            -0.6211
## ID4           -0.3521             0.34109            0.3411            -0.3521
## ID5           -0.4612            -0.46120            0.6374            -0.4612
##     Micrococcales Propionibacteriales Aeromonadales Alteromonadales 
## ID1         5.097               2.146       -0.5623           0.1309
## ID2         5.223               1.965        0.5790           1.0898
## ID3         1.864               1.576        2.7801           1.1706
## ID4         2.739               6.307       -0.3521          -0.3521
## ID5         2.104               1.148        1.6182           0.9251
##     Pasteurellales  Rhizobiales Rhodobacterales Oceanospirillales
## ID1         -0.5623      0.8240          2.9931           -0.5623
## ID2          0.1735      0.1735          0.5790            0.5790
## ID3         -0.6211      0.7652          1.8638            2.0179
## ID4         -0.3521     -0.3521          0.3411           -0.3521
## ID5         -0.4612      0.2320          0.9251            0.2320
##     Bifidobacteriales  Sphingomonadales Synergistales Tissierellales    NA_1
## ID1            -0.5623           0.1309       -0.5623          5.724 -0.5623
## ID2            -0.5197          -0.5197       -0.5197          2.119  0.1735
## ID3            -0.6211          -0.6211       -0.6211          8.315 -0.6211
## ID4            -0.3521          -0.3521       -0.3521          3.312 -0.3521
## ID5            -0.4612          -0.4612       -0.4612          1.331 -0.4612
##     Myxococcales Kineosporiales Frankiales Neisseriales  Xanthomonadales
## ID1       0.5364        0.13089    -0.5623       -0.5623          0.1309
## ID2       1.0898       -0.51966    -0.5197       -0.5197          0.1735
## ID3       0.9883        0.07201    -0.6211       -0.6211          3.7483
## ID4      -0.3521       -0.35205    -0.3521       -0.3521         -0.3521
## ID5      -0.4612        0.23195    -0.4612       -0.4612         -0.4612
##     Chromatiales Anaeroplasmatales  Caulobacterales Pseudomonadales_1
## ID1       0.8240           -0.56226         -0.5623           -0.5623
## ID2       1.2721           -0.51966          0.1735           -0.5197
## ID3       2.2121            0.07201         -0.6211           -0.6211
## ID4      -0.3521           -0.35205         -0.3521           -0.3521
## ID5       0.9251           -0.46120         -0.4612           -0.4612
##     Bacteriovoracales Verrucomicrobiales Bacteroidales Bdellovibrionales
## ID1           -0.5623            -0.5623       -0.5623            0.1309
## ID2           -0.5197            -0.5197       -0.5197           -0.5197
## ID3           -0.6211            -0.6211       -0.6211           -0.6211
## ID4           -0.3521            -0.3521       -0.3521            0.7466
## ID5           -0.4612            -0.4612       -0.4612           -0.4612
##     Enterobacterales  Flavobacteriales Erysipelotrichales Pirellulales
## ID1            1.2295         -0.56226            -0.5623      -0.5623
## ID2            0.5790          2.91433            -0.5197      -0.5197
## ID3            6.6036          0.07201            -0.6211      -0.6211
## ID4            0.7466         -0.35205            -0.3521      -0.3521
## ID5            6.7112         -0.46120            -0.4612      -0.4612
##     Holosporales Thermoanaerobacterales Caldisericales Campylobacterales
## ID1      -0.5623                 0.1309        -0.5623           -0.5623
## ID2      -0.5197                 0.1735        -0.5197           -0.5197
## ID3      -0.6211                -0.6211        -0.6211           -0.6211
## ID4      -0.3521                -0.3521        -0.3521           -0.3521
## ID5      -0.4612                -0.4612        -0.4612           -0.4612
##     Cardiobacteriales  Vibrionales  Catenulisporales      NA_2 CellVibrionales
## ID1           -0.56226       0.1309           -0.5623  0.13089        -0.56226
## ID2           -0.51966       0.1735           -0.5197  0.57895        -0.51966
## ID3            0.07201       1.1706           -0.6211  0.07201         0.07201
## ID4           -0.35205      -0.3521           -0.3521 -0.35205        -0.35205
## ID5           -0.46120      -0.4612           -0.4612 -0.46120         0.23195
##     Chitinispirillales  Chlorobiales Chroococcidiopsidales Cytophagales
## ID1             -0.5623      -0.5623               -0.5623      -0.5623
## ID2              0.1735      -0.5197               -0.5197      -0.5197
## ID3             -0.6211       0.4775               -0.6211      -0.6211
## ID4             -0.3521      -0.3521               -0.3521      -0.3521
## ID5             -0.4612      -0.4612               -0.4612      -0.4612
##     Coprothermobacterales Corynebacteriales Legionellales Parachlamydiales
## ID1               -0.5623             6.858       -0.5623          -0.5623
## ID2               -0.5197             5.177       -0.5197          -0.5197
## ID3               -0.6211             4.936       -0.6211          -0.6211
## ID4               -0.3521             6.326       -0.3521          -0.3521
## ID5               -0.4612             4.968       -0.4612          -0.4612
##     Cucurbitales Thiotrichales    NA_3 Dehalococcoidales Deinococcales 
## ID1      -0.5623        0.1309 -0.5623           -0.5623        -0.5623
## ID2      -0.5197       -0.5197 -0.5197           -0.5197        -0.5197
## ID3      -0.6211       -0.6211 -0.6211           -0.6211        -0.6211
## ID4      -0.3521       -0.3521 -0.3521            0.3411        -0.3521
## ID5      -0.4612       -0.4612 -0.4612           -0.4612        -0.4612
##     Burkholderiales  Syntrophobacterales DesulfoVibrionales Desulfobacterales 
## ID1          -0.5623              0.1309            -0.5623            -0.5623
## ID2           0.1735             -0.5197             0.1735            -0.5197
## ID3          -0.6211             -0.6211             0.7652            -0.6211
## ID4          -0.3521             -0.3521            -0.3521            -0.3521
## ID5          -0.4612             -0.4612            -0.4612            -0.4612
##     Desulfuromonadales Veillonellales Egibacterales Entomoplasmatales     NA_4
## ID1             0.1309        -0.5623       -0.5623           -0.5623 -0.56226
## ID2            -0.5197        -0.5197       -0.5197           -0.5197 -0.51966
## ID3            -0.6211        -0.6211       -0.6211           -0.6211  0.07201
## ID4             1.8452        -0.3521       -0.3521           -0.3521 -0.35205
## ID5            -0.4612        -0.4612       -0.4612            0.2320 -0.46120
##     Chitinophagales Cryptosporangiales Gaiellales  Gemmatimonadales
## ID1         -0.5623            -0.5623     -0.5623          -0.5623
## ID2         -0.5197            -0.5197      0.1735          -0.5197
## ID3         -0.6211            -0.6211     -0.6211          -0.6211
## ID4         -0.3521            -0.3521     -0.3521          -0.3521
## ID5         -0.4612            -0.4612     -0.4612          -0.4612
##     Geodermatophilales Micrococcales_1 Clostridiales  Rickettsiales
## ID1             0.5364         -0.5623        -0.5623       0.13089
## ID2            -0.5197         -0.5197        -0.5197      -0.51966
## ID3            -0.6211         -0.6211        -0.6211       0.07201
## ID4            -0.3521         -0.3521        -0.3521      -0.35205
## ID5            -0.4612         -0.4612        -0.4612      -0.46120
##     CellVibrionales  Glycomycetales Pasteurellales Halanaerobiales Jiangellales
## ID1          -0.5623        -0.5623        -0.5623         -0.5623      -0.5623
## ID2          -0.5197        -0.5197        -0.5197         -0.5197      -0.5197
## ID3          -0.6211        -0.6211        -0.6211         -0.6211      -0.6211
## ID4          -0.3521        -0.3521        -0.3521         -0.3521      -0.3521
## ID5          -0.4612        -0.4612        -0.4612         -0.4612      -0.4612
##     Synechococcales Haloplasmatales Holophagales    NA_5 Ignavibacteriales 
## ID1         -0.5623         -0.5623      -0.5623 -0.5623            -0.5623
## ID2         -0.5197         -0.5197      -0.5197 -0.5197            -0.5197
## ID3         -0.6211         -0.6211      -0.6211 -0.6211            -0.6211
## ID4         -0.3521         -0.3521      -0.3521 -0.3521            -0.3521
## ID5          0.2320         -0.4612      -0.4612 -0.4612            -0.4612
##     Kiloniellales Streptomycetales  Nautiliales  Limnochordales Mariprofundales
## ID1       -0.5623           -0.5623      -0.5623        -0.5623         -0.5623
## ID2       -0.5197           -0.5197      -0.5197        -0.5197         -0.5197
## ID3       -0.6211           -0.6211      -0.6211        -0.6211         -0.6211
## ID4       -0.3521           -0.3521      -0.3521        -0.3521         -0.3521
## ID5       -0.4612           -0.4612      -0.4612        -0.4612         -0.4612
##     Nitrosomonadales Methylococcales Mycoplasmatales Thermales Oligosphaerales
## ID1          0.13089         -0.5623         0.82403   -0.5623         -0.5623
## ID2         -0.51966          0.1735         0.57895   -0.5197         -0.5197
## ID3          0.07201         -0.6211         0.07201   -0.6211         -0.6211
## ID4         -0.35205         -0.3521        -0.35205   -0.3521         -0.3521
## ID5         -0.46120         -0.4612         0.92510   -0.4612         -0.4612
##     Coriobacteriales Anaerolineales Sphingobacteriales Puniceicoccales
## ID1          -0.5623        -0.5623             0.1309         -0.5623
## ID2          -0.5197        -0.5197            -0.5197         -0.5197
## ID3          -0.6211        -0.6211            -0.6211         -0.6211
## ID4          -0.3521        -0.3521            -0.3521         -0.3521
## ID5          -0.4612        -0.4612            -0.4612         -0.4612
##     Arenicellales Bacteroidia    NA_6 Eggerthellales Rubrobacterales
## ID1       -0.5623     -0.5623 -0.5623        -0.5623          0.1309
## ID2       -0.5197     -0.5197 -0.5197        -0.5197         -0.5197
## ID3       -0.6211     -0.6211 -0.6211        -0.6211         -0.6211
## ID4       -0.3521     -0.3521 -0.3521        -0.3521         -0.3521
## ID5       -0.4612     -0.4612 -0.4612        -0.4612         -0.4612
##     Sorangiineae Sphaeropleales Planctomycetales Spirochaetales Sporichthyales
## ID1      -0.5623        -0.5623          -0.5623        -0.5623        -0.5623
## ID2      -0.5197        -0.5197          -0.5197        -0.5197        -0.5197
## ID3      -0.6211         0.9883          -0.6211        -0.6211        -0.6211
## ID4      -0.3521        -0.3521          -0.3521        -0.3521        -0.3521
## ID5      -0.4612        -0.4612          -0.4612        -0.4612        -0.4612
##     Nevskiales Thermoanaerobaculales Nitrospirales Thermoflexales
## ID1    -0.5623               -0.5623       -0.5623        -0.5623
## ID2    -0.5197               -0.5197        0.1735        -0.5197
## ID3    -0.6211               -0.6211       -0.6211        -0.6211
## ID4    -0.3521               -0.3521       -0.3521        -0.3521
## ID5    -0.4612               -0.4612       -0.4612        -0.4612
##     Thermoleophilales Micrococcales_2 VampiroVibrionales Cystobacterineae diet
## ID1           -0.5623         -0.5623            -0.5623           0.1309  Veg
## ID2           -0.5197         -0.5197            -0.5197           0.5790  Veg
## ID3           -0.6211         -0.6211            -0.6211          -0.6211  Veg
## ID4           -0.3521         -0.3521            -0.3521          -0.3521  Veg
## ID5           -0.4612         -0.4612            -0.4612           0.6374  Veg
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
##      Mixed    14   9
##      Veg       9  15
##                                         
##                Accuracy : 0.617         
##                  95% CI : (0.464, 0.755)
##     No Information Rate : 0.511         
##     P-Value [Acc > NIR] : 0.0942        
##                                         
##                   Kappa : 0.234         
##                                         
##  Mcnemar's Test P-Value : 1.0000        
##                                         
##             Sensitivity : 0.609         
##             Specificity : 0.625         
##          Pos Pred Value : 0.609         
##          Neg Pred Value : 0.625         
##              Prevalence : 0.489         
##          Detection Rate : 0.298         
##    Detection Prevalence : 0.489         
##       Balanced Accuracy : 0.617         
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
tune_grid <- expand.grid(nrounds = c(100, 200),
                         max_depth = 6,
                         colsample_bytree = c(0.8, 1),
                         eta = c(0.1, 0.3),
                         gamma = 0,
                         min_child_weight = c(4, 5),
                         subsample = 0.8
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


```r
if( !require(plotROC) ){
    install.packages("plotROC")
    library(plotROC)
}
# Convert diets forms into binary
model$pred$pred_binary <- ifelse(model$pred$pred == "Mixed", 1, 0)
model$pred$obs_binary <- ifelse(model$pred$obs == "Mixed", 1, 0)

# Predictions
predictions <- model$pred

# Plot roc plot
plot_roc <- ggplot(predictions, aes_string(m = "pred_binary", d = "obs_binary")) + 
    geom_roc(labels = FALSE, n.cuts = 0)

# Create additional aesthetics and calculate AUC
plot <- plot_roc + 
    style_roc() + 
    geom_text(aes(0.5, 0.5, label = paste0("AUC = ", round(calc_auc(plot_roc)$AUC, 3) ) )) +
    theme(legend.position="none")
plot
```

![](40_machine_learning_files/figure-latex/super5-1.pdf)<!-- --> 


## Unsupervised machine learning

"Unsupervised" means that the labels (e.g., patient status is not known), 
and patterns are learned based only the abundance table, for instance. 
Unsupervised ML is also known as a data mining where patterns are extracted 
from big datasets. 

For unsupervised machine learning, please refer to chapters that are listed below:

- \@ref(clustering)
- \@ref(community-similarity) 

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
 [1] plotROC_2.3.0                  caret_6.0-93                  
 [3] lattice_0.20-45                ggplot2_3.3.6                 
 [5] mikropml_1.3.0                 mia_1.5.16                    
 [7] MultiAssayExperiment_1.22.0    TreeSummarizedExperiment_2.1.4
 [9] Biostrings_2.64.1              XVector_0.36.0                
[11] SingleCellExperiment_1.18.0    SummarizedExperiment_1.26.1   
[13] Biobase_2.56.0                 GenomicRanges_1.48.0          
[15] GenomeInfoDb_1.32.4            IRanges_2.30.1                
[17] S4Vectors_0.34.0               BiocGenerics_0.42.0           
[19] MatrixGenerics_1.8.1           matrixStats_0.62.0-9003       
[21] BiocStyle_2.24.0               rebook_1.6.0                  

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
 [29] crayon_1.5.2                RCurl_1.98-1.8             
 [31] jsonlite_1.8.0              graph_1.74.0               
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
 [61] XML_3.99-0.10               scuttle_1.6.3              
 [63] nnet_7.3-18                 CodeDepends_0.6.5          
 [65] utf8_1.2.2                  labeling_0.4.2             
 [67] tidyselect_1.1.2            rlang_1.0.6                
 [69] reshape2_1.4.4              munsell_0.5.0              
 [71] tools_4.2.1                 cachem_1.0.6               
 [73] xgboost_1.6.0.1             cli_3.4.1                  
 [75] DirichletMultinomial_1.38.0 generics_0.1.3             
 [77] RSQLite_2.2.17              evaluate_0.16              
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
[107] bitops_1.0-7                irlba_2.3.5                
[109] R6_2.5.1                    bookdown_0.29              
[111] gridExtra_2.3               vipor_0.4.5                
[113] parallelly_1.32.1           codetools_0.2-18           
[115] MASS_7.3-58.1               assertthat_0.2.1           
[117] withr_2.5.0                 GenomeInfoDbData_1.2.8     
[119] mgcv_1.8-40                 parallel_4.2.1             
[121] grid_4.2.1                  rpart_4.1.16               
[123] beachmat_2.12.0             timeDate_4021.104          
[125] tidyr_1.2.1                 class_7.3-20               
[127] rmarkdown_2.16              DelayedMatrixStats_1.18.1  
[129] pROC_1.18.0                 lubridate_1.8.0            
[131] ggbeeswarm_0.6.0           
```
</div>

