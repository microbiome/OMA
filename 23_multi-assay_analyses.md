# Multi-assay analyses {#multi-assay_analyses}

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

```
## Loading required package: ecodist
```


```r
library(mia)
```

Multi-omics means that we integrate data from multiple sources. For example, 
we can integrate microbial abundances in the gut with biomolecular profiling data
from blood samples. This kind of integrative multi-omic approaches can support the 
analysis of microbiome dysbiosis and facilitate the discovery of novel biomarkers 
for health and disease. 

With cross-correlation analysis, we can analyze how strongly and how differently
variables are associated between each other. For instance, we can analyze if 
higher presence of a specific taxon equals to higher levels of a biomolecule.

The data containers that the _miaverse_ utilizes are scalable and they can contain
different types of data in a same container. Because of that, the _miaverse_ is
well-suitable for multi-assay microbiome data which incorporates different types 
of complementary data sources in a single reproducible workflow. 

Another experiment can be stored in 
[altExp](https://microbiome.github.io/OMA/containers.html#alternative-experiments) 
slot of SE data container or both experiments can be stored side-by-side in
[MAE](https://microbiome.github.io/OMA/containers.html#multiassayexperiments) 
data container. 

Different experiments are first imported into SE or TreeSE data container similarly
to the case when only one experiment is present. After that different experiments are 
combined into the same data container. Result is one TreeSE object with alternative
experiment in altExp slot, or MAE object with multiple experiment in its 
experiment slot. 

As an example data, we use data from following publication: Hintikka L _et al._ (2021) 
[Xylo-oligosaccharides in prevention of hepatic steatosis and adipose tissue inflammation: 
associating taxonomic and metabolomic patterns in fecal microbiotas with 
biclustering](https://doi.org/10.3390/ijerph18084049).

In this article, mice were fed with high-fat and low-fat diets with or without prebiotics.
The purpose of this was to study if prebiotics would reduce the negative impacts
of high-fat diet. 

This example data can be loaded from microbiomeDataSets. The data is already in MAE
format. It includes three different experiments: microbial abundance data, 
metabolite concentrations, and data about different biomarkers. Help for importing
data into SE object you can find from [here](https://microbiome.github.io/OMA/containers.html#loading-experimental-microbiome-data).


```r
# Load the data
mae <- microbiomeDataSets::HintikkaXOData()

mae
```

```
## A MultiAssayExperiment object of 3 listed
##  experiments with user-defined names and respective classes.
##  Containing an ExperimentList class object of length 3:
##  [1] microbiota: SummarizedExperiment with 12706 rows and 40 columns
##  [2] metabolites: SummarizedExperiment with 38 rows and 40 columns
##  [3] biomarkers: SummarizedExperiment with 39 rows and 40 columns
## Functionality:
##  experiments() - obtain the ExperimentList instance
##  colData() - the primary/phenotype DataFrame
##  sampleMap() - the sample coordination DataFrame
##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
##  *Format() - convert into a long or wide DataFrame
##  assays() - convert ExperimentList to a SimpleList of matrices
##  exportClass() - save data to flat files
```

```r
if(!require(stringr)){
    install.packages("stringr")
    library(stringr)
}
# Drop off those bacteria that do not include information in Phylum or lower levels
mae[[1]] <- mae[[1]][!is.na(rowData(mae[[1]])$Phylum), ]
# Clean taxonomy data, so that names do not include additional characters
rowData(mae[[1]]) <- DataFrame(apply(rowData(mae[[1]]), 2, 
                                     str_remove, pattern = "._[0-9]__"))
# Microbiome data
mae[[1]]
```

```
## class: SummarizedExperiment 
## dim: 12613 40 
## metadata(0):
## assays(1): counts
## rownames(12613): GAYR01026362.62.2014 CVJT01000011.50.2173 ...
##   JRJTB:03787:02429 JRJTB:03787:02478
## rowData names(7): Phylum Class ... Species OTU
## colnames(40): C1 C2 ... C39 C40
## colData names(0):
```


```r
# Metabolite data
mae[[2]]
```

```
## class: SummarizedExperiment 
## dim: 38 40 
## metadata(0):
## assays(1): nmr
## rownames(38): Butyrate Acetate ... Malonate 1,3-dihydroxyacetone
## rowData names(0):
## colnames(40): C1 C2 ... C39 C40
## colData names(0):
```


```r
# Biomarker data
mae[[3]]
```

```
## class: SummarizedExperiment 
## dim: 39 40 
## metadata(0):
## assays(1): signals
## rownames(39): Triglycerides_liver CLSs_epi ... NPY_serum Glycogen_liver
## rowData names(0):
## colnames(40): C1 C2 ... C39 C40
## colData names(0):
```

## Cross-correlation Analysis

Next we can do the cross-correlation analysis. 
Here we analyse if individual bacteria genera correlates
with concentrations of individual metabolites. This helps as to answer the question: 
"If this bacteria is present, is this metabolite's concentration then low or high"?


```r
# Agglomerate microbiome data at Genus level
mae[[1]] <- agglomerateByRank(mae[[1]], rank = "Genus")
# Does log10 transform for microbiome data
mae[[1]] <- transformSamples(mae[[1]], method = "log10", pseudocount = 1)

# Cross correlates data sets
correlation_table <- testExperimentCrossCorrelation(mae, 
                                                    experiment1 = 1,
                                                    experiment2 = 2,
                                                    abund_values1 = "log10", 
                                                    abund_values2 = "nmr",
                                                    method = "spearman", 
                                                    p_adj_threshold = 0.05,
                                                    cor_threshold = 0)

knitr::kable(head(correlation_table))
```


\begin{tabular}{l|l|r|r|r}
\hline
Var1 & Var2 & cor & pval & p\_adj\\
\hline
Genus:Streptococcus & Butyrate & -0.4136 & 0.0080 & 0.1729\\
\hline
Genus:Ruminiclostridium 5 & Butyrate & 0.6760 & 0.0000 & 0.0017\\
\hline
Genus:Lactococcus & Butyrate & -0.3657 & 0.0209 & 0.2871\\
\hline
Genus:Globicatella & Butyrate & -0.4792 & 0.0018 & 0.0790\\
\hline
Genus:Enterobacter & Butyrate & -0.6011 & 0.0000 & 0.0079\\
\hline
Genus:Lachnoclostridium & Butyrate & -0.6346 & 0.0000 & 0.0039\\
\hline
\end{tabular}

Reorders the table


```r
if(!require(reshape2)){
    install.packages("reshape2")
}

# Converts data to matrix, correlations as values
mat <- reshape2::acast(correlation_table, Var2 ~ Var1, value.var = "cor")

# Hierarchical clustering, gets the order of taxa and lipids
taxa_indices <- hclust(as.dist(1 - cor(mat, use="pairwise.complete.obs")))$order
order_taxa <- colnames(mat)[taxa_indices]
lipids_indices <- hclust(as.dist(1 - cor(t(mat), use="pairwise.complete.obs")))$order
order_lipids <- rownames(mat)[lipids_indices]

# Converts taxa and lipids columns to factor so that they have the desired order
correlation_table[["Var1"]] <- factor(correlation_table[["Var1"]], levels = order_taxa)
correlation_table[["Var2"]] <- factor(correlation_table[["Var2"]], levels = order_lipids)
```

Creates the heatmap


```r
library(ggplot2)

# Determines the scaling of colors
limits <- c(-1, 1)
breaks <- seq(from = min(limits), to = max(limits), by = 0.2)
colours <- c("darkblue", "blue", "white", "red", "darkred")

# Which observation have p-value under 0.05? --> creates a subset
cor_table_sub <- correlation_table[which(correlation_table[["p_adj"]] < 0.05), ]

# Creates a ggplot object
ggplot(correlation_table, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile() +
  
  scale_fill_gradientn(name = "Correlation", 
                       breaks = breaks, limits = limits, colours = colours) + 
  
  # Adds label to those observations that have p-value under 0.05
  geom_text(data = cor_table_sub, aes(x = Var1, y = Var2, label = "+")) +
  
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.key.size = unit(1, "cm")) +
  labs(x = "Taxa", y = "Metabolites")
```

![](23_multi-assay_analyses_files/figure-latex/cross-correlation7-1.pdf)<!-- --> 

## Multi-Omics Factor Analysis

Multi-Omics Factor Analysis [@Argelaguet2018] (MOFA) is
an unsupervised method for integrating multi-omic data sets in a downstream analysis.
It could be seen as a generalization of principal component analysis. Yet, with the ability to infer a latent (low-dimensional) representation, 
shared among the mutliple (-omics) data sets in hand.

We use the R [MOFA2](https://biofam.github.io/MOFA2/index.html) package
for the analysis, and [install](https://biofam.github.io/MOFA2/installation.html) the corresponding dependencies.


```r
if(!require(MOFA2)){
    BiocManager::install("MOFA2")
}

# For inter-operability between Python and R, and setting Python dependencies,
# reticulate package is needed
if(!require(reticulate)){
    install.packages("reticulate")
}

reticulate::install_miniconda(force = TRUE)
```

```
## [1] "/github/home/.local/share/r-miniconda"
```

```r
reticulate::use_miniconda(condaenv = "env1", required = FALSE)
reticulate::py_install(packages = c("mofapy2"), pip = TRUE)
```

The `mae` object could be used straight to create the MOFA model. Yet, we transform 
our assays since the model assumes normality per default. Other distributions that
can be used, include Poisson or Bernoulli.


```r
library(MOFA2)
# For simplicity, classify all high-fat diets as high-fat, and all the low-fat 
# diets as low-fat diets
colData(mae)$Diet <- ifelse(colData(mae)$Diet == "High-fat" | 
                              colData(mae)$Diet == "High-fat + XOS", 
                            "High-fat", "Low-fat")

# Removing duplicates at the microbiome data
# which are also in form e.g. "Ambiguous" and "uncultured" taxa
mae[[1]] <- mae[[1]][!duplicated(rownames(assay(mae[[1]]))), ]

# Transforming microbiome data with rclr
mae[[1]] <- transformCounts(mae[[1]], abund_values = "counts", method = "rclr")

# Transforming metabolomic data with log10
mae[[2]] <- transformSamples(mae[[2]], abund_values = "nmr", method = "log10")

# Transforming biomarker data with z-transform
mae[[3]] <- transformFeatures(mae[[3]], abund_values = "signals", method = "z", pseudocount = 1)

# Removing assays no longer needed
assay(mae[[1]], "counts") <- NULL
assay(mae[[1]], "log10") <- NULL
assay(mae[[2]], "nmr") <- NULL
assay(mae[[3]], "signals") <- NULL

# Building our mofa model
model <- create_mofa_from_MultiAssayExperiment(mae,
                                               groups = "Diet", 
                                               extract_metadata = TRUE)
model
```

```
## Untrained MOFA model with the following characteristics: 
##  Number of views: 3 
##  Views names: microbiota metabolites biomarkers 
##  Number of features (per view): 239 38 39 
##  Number of groups: 2 
##  Groups names: High-fat Low-fat 
##  Number of samples (per group): 20 20 
## 
```

Model options could be defined as follows:


```r
model_opts <- get_default_model_options(model)
model_opts$num_factors <- 15
head(model_opts)
```

```
## $likelihoods
##  microbiota metabolites  biomarkers 
##  "gaussian"  "gaussian"  "gaussian" 
## 
## $num_factors
## [1] 15
## 
## $spikeslab_factors
## [1] FALSE
## 
## $spikeslab_weights
## [1] TRUE
## 
## $ard_factors
## [1] TRUE
## 
## $ard_weights
## [1] TRUE
```

Model's training options are defined with the following:


```r
train_opts <- get_default_training_options(model)
head(train_opts)
```

```
## $maxiter
## [1] 1000
## 
## $convergence_mode
## [1] "fast"
## 
## $drop_factor_threshold
## [1] -1
## 
## $verbose
## [1] FALSE
## 
## $startELBO
## [1] 1
## 
## $freqELBO
## [1] 5
```

Preparing and training the model:


```r
model.prepared <- prepare_mofa(
  object = model,
  model_options = model_opts
)
model.trained <- run_mofa(model.prepared)
```

Visualizing the variance explained:


```r
library(patchwork)
wrap_plots(
    plot_variance_explained(model.trained, x="view", y="factor", plot_total = T),
    nrow = 2
) + plot_annotation(title = "Variance Explained per factor and assay",
                    theme = theme(plot.title = element_text(hjust = 0.5)))
```

![](23_multi-assay_analyses_files/figure-latex/unnamed-chunk-6-1.pdf)<!-- --> 

The top weights for each assay using the three first factors:


```r
plots <- lapply(c("microbiota", "metabolites","biomarkers"), function(name) {
    plot_top_weights(model.trained,
                     view = name,
                     factors = 1:3,
                     nfeatures = 10) +
        labs(title = paste0("Top weights of the ", name," assay"))
})
wrap_plots(plots, nrow = 3) & theme(text = element_text(size = 8))
```

![](23_multi-assay_analyses_files/figure-latex/unnamed-chunk-7-1.pdf)<!-- --> 

More tutorials and examples of using the package are found at: [link](https://biofam.github.io/MOFA2/tutorials.html)



## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

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
 [1] patchwork_1.1.1                reticulate_1.24               
 [3] MOFA2_1.4.0                    ggplot2_3.3.5                 
 [5] reshape2_1.4.4                 stringr_1.4.0                 
 [7] microbiomeDataSets_1.1.5       mia_1.3.14                    
 [9] MultiAssayExperiment_1.20.0    TreeSummarizedExperiment_2.1.4
[11] Biostrings_2.62.0              XVector_0.34.0                
[13] SingleCellExperiment_1.16.0    SummarizedExperiment_1.24.0   
[15] Biobase_2.54.0                 GenomicRanges_1.46.1          
[17] GenomeInfoDb_1.30.1            IRanges_2.28.0                
[19] S4Vectors_0.32.3               BiocGenerics_0.40.0           
[21] MatrixGenerics_1.6.0           matrixStats_0.61.0-9001       
[23] ecodist_2.0.7                  BiocStyle_2.22.0              
[25] rebook_1.4.0                  

loaded via a namespace (and not attached):
  [1] AnnotationHub_3.2.1           corrplot_0.92                
  [3] BiocFileCache_2.2.1           plyr_1.8.6                   
  [5] lazyeval_0.2.2                splines_4.1.2                
  [7] BiocParallel_1.28.3           scater_1.22.0                
  [9] digest_0.6.29                 yulab.utils_0.0.4            
 [11] htmltools_0.5.2               viridis_0.6.2                
 [13] fansi_1.0.2                   magrittr_2.0.2               
 [15] memoise_2.0.1                 ScaledMatrix_1.2.0           
 [17] cluster_2.1.2                 DECIPHER_2.22.0              
 [19] colorspace_2.0-2              blob_1.2.2                   
 [21] rappdirs_0.3.3                ggrepel_0.9.1                
 [23] xfun_0.29                     dplyr_1.0.7                  
 [25] crayon_1.4.2                  RCurl_1.98-1.5               
 [27] jsonlite_1.7.3                graph_1.72.0                 
 [29] ape_5.6-1                     glue_1.6.1                   
 [31] gtable_0.3.0                  zlibbioc_1.40.0              
 [33] DelayedArray_0.20.0           BiocSingular_1.10.0          
 [35] Rhdf5lib_1.16.0               HDF5Array_1.22.1             
 [37] scales_1.1.1                  pheatmap_1.0.12              
 [39] DBI_1.1.2                     Rcpp_1.0.8                   
 [41] viridisLite_0.4.0             xtable_1.8-4                 
 [43] decontam_1.14.0               tidytree_0.3.7               
 [45] bit_4.0.4                     rsvd_1.0.5                   
 [47] httr_1.4.2                    RColorBrewer_1.1-2           
 [49] dir.expiry_1.2.0              ellipsis_0.3.2               
 [51] farver_2.1.0                  pkgconfig_2.0.3              
 [53] XML_3.99-0.8                  scuttle_1.4.0                
 [55] uwot_0.1.11                   CodeDepends_0.6.5            
 [57] dbplyr_2.1.1                  here_1.0.1                   
 [59] utf8_1.2.2                    labeling_0.4.2               
 [61] tidyselect_1.1.1              rlang_1.0.0                  
 [63] later_1.3.0                   AnnotationDbi_1.56.2         
 [65] munsell_0.5.0                 BiocVersion_3.14.0           
 [67] tools_4.1.2                   cachem_1.0.6                 
 [69] cli_3.1.1                     DirichletMultinomial_1.36.0  
 [71] generics_0.1.2                RSQLite_2.2.9                
 [73] ExperimentHub_2.2.1           evaluate_0.14                
 [75] fastmap_1.1.0                 yaml_2.2.2                   
 [77] knitr_1.37                    bit64_4.0.5                  
 [79] purrr_0.3.4                   KEGGREST_1.34.0              
 [81] nlme_3.1-155                  sparseMatrixStats_1.6.0      
 [83] mime_0.12                     compiler_4.1.2               
 [85] beeswarm_0.4.0                filelock_1.0.2               
 [87] curl_4.3.2                    png_0.1-7                    
 [89] interactiveDisplayBase_1.32.0 treeio_1.18.1                
 [91] tibble_3.1.6                  stringi_1.7.6                
 [93] basilisk.utils_1.6.0          highr_0.9                    
 [95] forcats_0.5.1                 lattice_0.20-45              
 [97] Matrix_1.4-0                  vegan_2.5-7                  
 [99] permute_0.9-7                 vctrs_0.3.8                  
[101] rhdf5filters_1.6.0            pillar_1.7.0                 
[103] lifecycle_1.0.1               BiocManager_1.30.16          
[105] BiocNeighbors_1.12.0          cowplot_1.1.1                
[107] bitops_1.0-7                  irlba_2.3.5                  
[109] httpuv_1.6.5                  R6_2.5.1                     
[111] bookdown_0.24                 promises_1.2.0.1             
[113] gridExtra_2.3                 vipor_0.4.5                  
[115] codetools_0.2-18              MASS_7.3-55                  
[117] assertthat_0.2.1              rhdf5_2.38.0                 
[119] rprojroot_2.0.2               withr_2.4.3                  
[121] GenomeInfoDbData_1.2.7        mgcv_1.8-38                  
[123] parallel_4.1.2                grid_4.1.2                   
[125] beachmat_2.10.0               basilisk_1.6.0               
[127] tidyr_1.2.0                   rmarkdown_2.11               
[129] DelayedMatrixStats_1.16.0     Rtsne_0.15                   
[131] shiny_1.7.1                   ggbeeswarm_0.6.0             
```
</div>

