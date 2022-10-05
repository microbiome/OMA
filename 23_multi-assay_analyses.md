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
# Agglomerate microbiome data at family level
mae[[1]] <- agglomerateByPrevalence(mae[[1]], rank = "Family")
# Does log10 transform for microbiome data
mae[[1]] <- transformSamples(mae[[1]], method = "log10", pseudocount = 1)

# Give unique names so that we do not have problems when we are creating a plot
rownames(mae[[1]]) <- getTaxonomyLabels(mae[[1]])

# Cross correlates data sets
correlations <- testExperimentCrossCorrelation(mae, 
                                               experiment1 = 1,
                                               experiment2 = 2,
                                               assay_name1 = "log10", 
                                               assay_name2 = "nmr",
                                               method = "spearman", 
                                               p_adj_threshold = NULL,
                                               cor_threshold = NULL,
                                               # Remove when mia is fixed
                                               mode = "matrix",
                                               sort = TRUE,
                                               show_warnings = FALSE)
```

Creates the heatmap


```r
if( !require("ComplexHeatmap") ){
    BiocManager::install("ComplexHeatmap")
    library("ComplexHeatmap")
}

# Create a heatmap and store it
plot <- Heatmap(correlations$cor,
                # Print values to cells
                cell_fun = function(j, i, x, y, width, height, fill) {
                    # If the p-value is under threshold
                    if( !is.na(correlations$p_adj[i, j]) & correlations$p_adj[i, j] < 0.05 ){
                        # Print "X"
                        grid.text(sprintf("%s", "X"), x, y, gp = gpar(fontsize = 8, col = "black"))
                        }
                    },
                heatmap_legend_param = list(title = "", legend_height = unit(5, "cm"))
                )
plot
```

![](23_multi-assay_analyses_files/figure-latex/cross-correlation6-1.pdf)<!-- --> 

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
mae[[1]] <- transformCounts(mae[[1]], method = "relabundance")
mae[[1]] <- transformCounts(mae[[1]], assay_name = "relabundance", method = "rclr")

# Transforming metabolomic data with log10
mae[[2]] <- transformSamples(mae[[2]], assay_name = "nmr", method = "log10")

# Transforming biomarker data with z-transform
mae[[3]] <- transformFeatures(mae[[3]], assay_name = "signals", method = "z", pseudocount = 1)

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
##  Number of features (per view): 38 38 39 
##  Number of groups: 2 
##  Groups names: High-fat Low-fat 
##  Number of samples (per group): 20 20 
## 
```

Model options could be defined as follows:


```r
model_opts <- get_default_model_options(model)
model_opts$num_factors <- 5
head(model_opts)
```

```
## $likelihoods
##  microbiota metabolites  biomarkers 
##  "gaussian"  "gaussian"  "gaussian" 
## 
## $num_factors
## [1] 5
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
library(ggplot2)
wrap_plots(
    plot_variance_explained(model.trained, x="view", y="factor", plot_total = T),
    nrow = 2
) + plot_annotation(title = "Variance Explained per factor and assay",
                    theme = theme(plot.title = element_text(hjust = 0.5)))
```

![](23_multi-assay_analyses_files/figure-latex/unnamed-chunk-6-1.pdf)<!-- --> 

The top weights for each assay using all 5 factors:


```r
plots <- lapply(c("microbiota", "metabolites","biomarkers"), function(name) {
    plot_top_weights(model.trained,
                     view = name,
                     factors = "all",
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
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] ggplot2_3.3.6                  patchwork_1.1.2               
 [3] reticulate_1.26                MOFA2_1.6.0                   
 [5] ComplexHeatmap_2.12.1          stringr_1.4.1                 
 [7] microbiomeDataSets_1.1.7       mia_1.5.16                    
 [9] MultiAssayExperiment_1.22.0    TreeSummarizedExperiment_2.1.4
[11] Biostrings_2.64.1              XVector_0.36.0                
[13] SingleCellExperiment_1.18.1    SummarizedExperiment_1.26.1   
[15] Biobase_2.56.0                 GenomicRanges_1.48.0          
[17] GenomeInfoDb_1.32.4            IRanges_2.30.1                
[19] S4Vectors_0.34.0               BiocGenerics_0.42.0           
[21] MatrixGenerics_1.8.1           matrixStats_0.62.0-9003       
[23] BiocStyle_2.24.0               rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] circlize_0.4.15               AnnotationHub_3.4.0          
  [3] corrplot_0.92                 BiocFileCache_2.4.0          
  [5] plyr_1.8.7                    lazyeval_0.2.2               
  [7] splines_4.2.1                 BiocParallel_1.30.3          
  [9] scater_1.24.0                 digest_0.6.29                
 [11] foreach_1.5.2                 yulab.utils_0.0.5            
 [13] htmltools_0.5.3               viridis_0.6.2                
 [15] fansi_1.0.3                   magrittr_2.0.3               
 [17] memoise_2.0.1                 ScaledMatrix_1.4.1           
 [19] cluster_2.1.4                 doParallel_1.0.17            
 [21] DECIPHER_2.24.0               colorspace_2.0-3             
 [23] blob_1.2.3                    rappdirs_0.3.3               
 [25] ggrepel_0.9.1                 xfun_0.33                    
 [27] dplyr_1.0.10                  crayon_1.5.2                 
 [29] RCurl_1.98-1.9                jsonlite_1.8.2               
 [31] graph_1.74.0                  iterators_1.0.14             
 [33] ape_5.6-2                     glue_1.6.2                   
 [35] gtable_0.3.1                  zlibbioc_1.42.0              
 [37] GetoptLong_1.0.5              DelayedArray_0.22.0          
 [39] BiocSingular_1.12.0           Rhdf5lib_1.18.2              
 [41] shape_1.4.6                   HDF5Array_1.24.2             
 [43] scales_1.2.1                  pheatmap_1.0.12              
 [45] DBI_1.1.3                     Rcpp_1.0.9                   
 [47] viridisLite_0.4.1             xtable_1.8-4                 
 [49] clue_0.3-61                   decontam_1.16.0              
 [51] tidytree_0.4.1                bit_4.0.4                    
 [53] rsvd_1.0.5                    httr_1.4.4                   
 [55] RColorBrewer_1.1-3            dir.expiry_1.4.0             
 [57] ellipsis_0.3.2                farver_2.1.1                 
 [59] pkgconfig_2.0.3               XML_3.99-0.11                
 [61] scuttle_1.6.3                 uwot_0.1.14                  
 [63] CodeDepends_0.6.5             dbplyr_2.2.1                 
 [65] here_1.0.1                    utf8_1.2.2                   
 [67] labeling_0.4.2                tidyselect_1.1.2             
 [69] rlang_1.0.6                   reshape2_1.4.4               
 [71] later_1.3.0                   AnnotationDbi_1.58.0         
 [73] munsell_0.5.0                 BiocVersion_3.15.2           
 [75] tools_4.2.1                   cachem_1.0.6                 
 [77] cli_3.4.1                     DirichletMultinomial_1.38.0  
 [79] generics_0.1.3                RSQLite_2.2.17               
 [81] ExperimentHub_2.4.0           evaluate_0.16                
 [83] fastmap_1.1.0                 yaml_2.3.5                   
 [85] knitr_1.40                    bit64_4.0.5                  
 [87] purrr_0.3.4                   KEGGREST_1.36.3              
 [89] nlme_3.1-159                  sparseMatrixStats_1.8.0      
 [91] mime_0.12                     compiler_4.2.1               
 [93] beeswarm_0.4.0                filelock_1.0.2               
 [95] curl_4.3.2                    png_0.1-7                    
 [97] interactiveDisplayBase_1.34.0 treeio_1.20.2                
 [99] tibble_3.1.8                  stringi_1.7.8                
[101] basilisk.utils_1.8.0          highr_0.9                    
[103] forcats_0.5.2                 lattice_0.20-45              
[105] Matrix_1.5-1                  vegan_2.6-2                  
[107] permute_0.9-7                 vctrs_0.4.2                  
[109] rhdf5filters_1.8.0            pillar_1.8.1                 
[111] lifecycle_1.0.2               BiocManager_1.30.18          
[113] GlobalOptions_0.1.2           BiocNeighbors_1.14.0         
[115] cowplot_1.1.1                 bitops_1.0-7                 
[117] irlba_2.3.5.1                 httpuv_1.6.6                 
[119] R6_2.5.1                      bookdown_0.29                
[121] promises_1.2.0.1              gridExtra_2.3                
[123] vipor_0.4.5                   codetools_0.2-18             
[125] MASS_7.3-58.1                 assertthat_0.2.1             
[127] rhdf5_2.40.0                  rprojroot_2.0.3              
[129] rjson_0.2.21                  withr_2.5.0                  
[131] GenomeInfoDbData_1.2.8        mgcv_1.8-40                  
[133] parallel_4.2.1                beachmat_2.12.0              
[135] basilisk_1.8.1                tidyr_1.2.1                  
[137] rmarkdown_2.16                DelayedMatrixStats_1.18.1    
[139] Rtsne_0.16                    Cairo_1.6-0                  
[141] shiny_1.7.2                   ggbeeswarm_0.6.0             
```
</div>

