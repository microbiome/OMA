# Differential Abundance {#differential-abundance}

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



A number of methods for differential abundance analysis are available,
and reviewed elsewhere.

## Tree-based methods

### Group-wise associations testing based on balances

[TreeSummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html) 
frequently includes a Phylogenetic tree along with associated data about the experiment (at `colData`), that holds covariates which can be used for analyzing group-wise associations. 

Such an analysis could be performed with the function `pibble` from the `fido` package, that offers a Multinomial Logistic-Normal Linear Regression model; see [vignette](https://jsilve24.github.io/fido/articles/introduction-to-fido.html) of package.

The following presents such an exemplary analysis based on the [Sprockett et al. (2020)](https://doi.org/10.1038/s41467-020-17541-6) available through `microbiomeDataSets` package.


```r
if (!require(fido)){
  # installing the fido package
  devtools::install_github("jsilve24/fido")
}
```

```
## driver       (NA -> 16e449946...) [GitHub]
## distribut... (NA -> 0.2.2       ) [CRAN]
## tensorA      (NA -> 0.36.2      ) [CRAN]
## checkmate    (NA -> 2.0.0       ) [CRAN]
## svUnit       (NA -> 1.0.6       ) [CRAN]
## HDInterval   (NA -> 0.2.2       ) [CRAN]
## RcppGSL      (NA -> 0.3.9       ) [CRAN]
## posterior    (NA -> 1.0.1       ) [CRAN]
## arrayhelpers (NA -> 1.1-0       ) [CRAN]
## coda         (NA -> 0.19-4      ) [CRAN]
## ggdist       (NA -> 3.0.0       ) [CRAN]
## RcppZiggurat (NA -> 0.1.6       ) [CRAN]
## RcppNumer... (NA -> 0.4-0       ) [CRAN]
## tidybayes    (NA -> 3.0.0       ) [CRAN]
## 
##      checking for file ‘/tmp/RtmpiALWgt/remotes7c329de8600/jsilve24-driver-16e4499/DESCRIPTION’ ...  ✔  checking for file ‘/tmp/RtmpiALWgt/remotes7c329de8600/jsilve24-driver-16e4499/DESCRIPTION’
##   ─  preparing ‘driver’:
##   ✔  checking DESCRIPTION meta-information
##   ─  checking for LF line-endings in source and make files and shell scripts
##   ─  checking for empty or unneeded directories
##      Omitted ‘LazyData’ from DESCRIPTION
##   ─  building ‘driver_0.1.1.tar.gz’
##      
##      checking for file ‘/tmp/RtmpiALWgt/remotes7c332b730a4/jsilve24-fido-c692141/DESCRIPTION’ ...  ✔  checking for file ‘/tmp/RtmpiALWgt/remotes7c332b730a4/jsilve24-fido-c692141/DESCRIPTION’
##   ─  preparing ‘fido’:
##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
##   ─  cleaning src
## ─  running ‘cleanup’
##      Warning: /tmp/RtmpdNtyvB/Rbuild95fe7f9cc4/fido/man/loglikPibbleCollapsed.Rd:28: unknown macro '\item'
##    Warning: /tmp/RtmpdNtyvB/Rbuild95fe7f9cc4/fido/man/loglikPibbleCollapsed.Rd:30: unknown macro '\item'
##    Warning: /tmp/RtmpdNtyvB/Rbuild95fe7f9cc4/fido/man/loglikPibbleCollapsed.Rd:33: unexpected section header '\value'
##      Warning: /tmp/RtmpdNtyvB/Rbuild95fe7f9cc4/fido/man/loglikPibbleCollapsed.Rd:41: unexpected section header '\description'
##    Warning: /tmp/RtmpdNtyvB/Rbuild95fe7f9cc4/fido/man/loglikPibbleCollapsed.Rd:51: unexpected section header '\examples'
##    Warning: /tmp/RtmpdNtyvB/Rbuild95fe7f9cc4/fido/man/loglikPibbleCollapsed.Rd:83: unexpected END_OF_INPUT '
##    '
##   ─  checking for LF line-endings in source and make files and shell scripts
##   ─  checking for empty or unneeded directories
##   ─  building ‘fido_0.1.13.tar.gz’
##      
## 
```

Loading the libraries and importing data:


```r
library(fido)
library(mia)
library(microbiomeDataSets)

tse <- SprockettTHData()
names(colData(tse))
```

```
##  [1] "Sample_ID"                "Sample_Type"             
##  [3] "Subject_ID"               "Sex"                     
##  [5] "Country"                  "Age_Days"                
##  [7] "Age_Months"               "Age_Years"               
##  [9] "Days_After_Birth"         "Months_After_Birth"      
## [11] "Age_Group"                "Age_Class"               
## [13] "Feeding_Status_Corrected" "Study_Code"              
## [15] "AgeGroup_By_SampleType"   "Dyad"                    
## [17] "Delivery_Mode"            "Community"               
## [19] "Community_Type"           "Community_Subtype"       
## [21] "Sample_Label"             "Cohort"                  
## [23] "Sample_Collection_Date"   "Community_Alias"
```

We pick three covariates ("Sex","Age_Years","Delivery_Mode") during this analysis as an example, and beforehand we check for missing data:


```r
cov_names <- c("Sex","Age_Years","Delivery_Mode")
na_counts <- apply(is.na(colData(tse)[,cov_names]), 2, sum)
na_summary<-as.data.frame(na_counts,row.names=cov_names)
na_summary
```

```
##               na_counts
## Sex                   0
## Age_Years          1603
## Delivery_Mode      1100
```

We drop samples with na values at the covariates (features) under analysis: 


```r
tse <- tse[ , !is.na(colData(tse)$Delivery_Mode) ]
tse <- tse[ , !is.na(colData(tse)$Age_Years) ]
tse
```

```
## class: TreeSummarizedExperiment 
## dim: 2319 483 
## metadata(0):
## assays(1): counts
## rownames(2319): ASV_4_Candidatus_Methanogranum ASV_9_Methanobrevibacter
##   ... ASV_6837_Proteus_cibarius/hauseri+4 ASV_6838_Proteus
## rowData names(8): Kingdom Phylum ... Species lowest_rank
## colnames(483): 100.BF1 101.BF9A ... 98.MF35A 99.MF134A
## colData names(24): Sample_ID Sample_Type ... Sample_Collection_Date
##   Community_Alias
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (2319 rows)
## rowTree: 1 phylo tree(s) (2319 leaves)
## colLinks: NULL
## colTree: NULL
## referenceSeq: a DNAStringSet (2319 sequences)
```

We agglomerate the data at a Phylum rank.
Note: Large assay data (along with the covariates/features data) could prevent the analysis later,
since the computation will construct matrices that would not always fit memory.


```r
tse_phylum <- agglomerateByRank(tse, "Phylum")
tse_phylum
```

```
## class: TreeSummarizedExperiment 
## dim: 20 483 
## metadata(0):
## assays(1): counts
## rownames(20): Phylum:Euryarchaeota Phylum:Firmicutes ...
##   Phylum:Fusobacteria Phylum:Epsilonbacteraeota
## rowData names(8): Kingdom Phylum ... Species lowest_rank
## colnames(483): 100.BF1 101.BF9A ... 98.MF35A 99.MF134A
## colData names(24): Sample_ID Sample_Type ... Sample_Collection_Date
##   Community_Alias
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (20 rows)
## rowTree: 1 phylo tree(s) (2319 leaves)
## colLinks: NULL
## colTree: NULL
## referenceSeq: a DNAStringSet (20 sequences)
```

We extract the counts assay and feature data to build the model matrix having an extra row of ones presenting the intercept for the regression task later: 


```r
Y <- assays(tse_phylum)$counts
# design matrix
# taking 3 covariates
sample_data<-as.data.frame(colData(tse_phylum)[,cov_names])
X <- t(model.matrix(~Sex+Age_Years+Delivery_Mode,data=sample_data))
X[,1:5]
```

```
##                      100.BF1 101.BF9A 102.BF36A 103.BF53A 104.BF83A
## (Intercept)             1.00     1.00      1.00      1.00      1.00
## SexMale                 0.00     0.00      0.00      0.00      0.00
## Age_Years               0.88     0.89      0.94      0.99      1.05
## Delivery_ModeVaginal    1.00     1.00      1.00      1.00      1.00
```

```r
Y[1:5,1:5]
```

```
##                       100.BF1 101.BF9A 102.BF36A 103.BF53A 104.BF83A
## Phylum:Euryarchaeota        0        0         0         0         0
## Phylum:Firmicutes       10724    12129     15453     14446      4160
## Phylum:Tenericutes          0      173         0         3         2
## Phylum:Actinobacteria   18127    19122      1395      9666     19633
## Kingdom:Eukaryota           0        0         0         0         0
```

Building the parameters for the `pibble` call to build the model; see more at [vignette](https://jsilve24.github.io/fido/articles/introduction-to-fido.html):


```r
n_taxa<-nrow(Y)
upsilon <- n_taxa+3
Omega <- diag(n_taxa)
G <- cbind(diag(n_taxa-1), -1)
Xi <- (upsilon-n_taxa)*G%*%Omega%*%t(G)
Theta <- matrix(0, n_taxa-1, nrow(X))
Gamma <- diag(nrow(X))
```

Automatically initializing the priors and visualizing their distributions:


```r
priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)
names_covariates(priors) <- rownames(X)
fido::plot(priors, pars="Lambda") + ggplot2::xlim(c(-5, 5))
```

<img src="17_differential_abundance_files/figure-html/unnamed-chunk-8-1.png" width="672" />

Estimating the posterior by including the data at `Y`.
Note: Some computational failures could occur (see [discussion](https://github-wiki-see.page/m/jsilve24/fido/wiki/Frequently-Asked-Questions))
the arguments `multDirichletBoot` `calcGradHess` could be passed in such case.


```r
priors$Y <- Y 
posterior <- refit(priors, optim_method="adam", multDirichletBoot=0.5) # ,, calcGradHess=FALSE
```

Printing a summary about the posterior predictive distribution:


```r
ppc_summary(posterior)
```

```
## Proportions of Observations within 95% Credible Interval: 0.9981
```
Plotting the summary of the posterior distributions of the regression parameters:


```r
names_categories(posterior) <- rownames(Y)
fido::plot(posterior,par="Lambda",focus.cov=rownames(X)[2:4])
```

<img src="17_differential_abundance_files/figure-html/unnamed-chunk-11-1.png" width="672" />

Seemingly the covariate "Age_Years" does not have effect on the model as "Delivery_Mode" would,
and "Sex" to some extent. Let's take a closer look at the two latter ones:


```r
fido::plot(posterior, par="Lambda", focus.cov = rownames(X)[c(2,4)])
```

<img src="17_differential_abundance_files/figure-html/unnamed-chunk-12-1.png" width="672" />

## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] microbiomeDataSets_1.1.1       MultiAssayExperiment_1.19.5   
 [3] mia_1.1.8                      TreeSummarizedExperiment_2.1.3
 [5] Biostrings_2.61.2              XVector_0.33.0                
 [7] SingleCellExperiment_1.15.1    SummarizedExperiment_1.23.1   
 [9] Biobase_2.53.0                 GenomicRanges_1.45.0          
[11] GenomeInfoDb_1.29.3            IRanges_2.27.0                
[13] S4Vectors_0.31.0               BiocGenerics_0.39.1           
[15] MatrixGenerics_1.5.2           matrixStats_0.60.0            
[17] fido_0.1.13                    BiocStyle_2.21.3              
[19] rebook_1.3.0                  

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                    tidyselect_1.1.1             
  [3] RSQLite_2.2.7                 AnnotationDbi_1.55.1         
  [5] grid_4.1.0                    BiocParallel_1.27.2          
  [7] devtools_2.4.2                munsell_0.5.0                
  [9] ScaledMatrix_1.1.0            codetools_0.2-18             
 [11] withr_2.4.2                   colorspace_2.0-2             
 [13] filelock_1.0.2                highr_0.9                    
 [15] knitr_1.33                    labeling_0.4.2               
 [17] GenomeInfoDbData_1.2.6        bit64_4.0.5                  
 [19] farver_2.1.0                  driver_0.1.1                 
 [21] rprojroot_2.0.2               coda_0.19-4                  
 [23] vctrs_0.3.8                   treeio_1.17.2                
 [25] generics_0.1.0                xfun_0.24                    
 [27] BiocFileCache_2.1.1           R6_2.5.0                     
 [29] ggbeeswarm_0.6.0              rsvd_1.0.5                   
 [31] bitops_1.0-7                  cachem_1.0.5                 
 [33] DelayedArray_0.19.1           assertthat_0.2.1             
 [35] promises_1.2.0.1              scales_1.1.1                 
 [37] beeswarm_0.4.0                gtable_0.3.0                 
 [39] beachmat_2.9.0                processx_3.5.2               
 [41] rlang_0.4.11                  splines_4.1.0                
 [43] lazyeval_0.2.2                checkmate_2.0.0              
 [45] BiocManager_1.30.16           yaml_2.2.1                   
 [47] reshape2_1.4.4                abind_1.4-5                  
 [49] backports_1.2.1               httpuv_1.6.1                 
 [51] tensorA_0.36.2                tools_4.1.0                  
 [53] usethis_2.0.1                 bookdown_0.22                
 [55] ggplot2_3.3.5                 ellipsis_0.3.2               
 [57] decontam_1.13.0               RColorBrewer_1.1-2           
 [59] jquerylib_0.1.4               posterior_1.0.1              
 [61] sessioninfo_1.1.1             Rcpp_1.0.7                   
 [63] plyr_1.8.6                    sparseMatrixStats_1.5.0      
 [65] zlibbioc_1.39.0               purrr_0.3.4                  
 [67] RCurl_1.98-1.3                ps_1.6.0                     
 [69] prettyunits_1.1.1             viridis_0.6.1                
 [71] ggrepel_0.9.1                 cluster_2.1.2                
 [73] fs_1.5.0                      DECIPHER_2.21.0              
 [75] magrittr_2.0.1                ggdist_3.0.0                 
 [77] tidybayes_3.0.0               pkgload_1.2.1                
 [79] mime_0.11                     evaluate_0.14                
 [81] arrayhelpers_1.1-0            xtable_1.8-4                 
 [83] XML_3.99-0.6                  gridExtra_2.3                
 [85] testthat_3.0.4                compiler_4.1.0               
 [87] scater_1.21.3                 tibble_3.1.3                 
 [89] crayon_1.4.1                  htmltools_0.5.1.1            
 [91] mgcv_1.8-36                   later_1.2.0                  
 [93] tidyr_1.1.3                   DBI_1.1.1                    
 [95] ExperimentHub_2.1.4           dbplyr_2.1.1                 
 [97] MASS_7.3-54                   rappdirs_0.3.3               
 [99] Matrix_1.3-4                  permute_0.9-5                
[101] cli_3.0.1                     parallel_4.1.0               
[103] pkgconfig_2.0.3               dir.expiry_1.1.0             
[105] scuttle_1.3.0                 svUnit_1.0.6                 
[107] vipor_0.4.5                   bslib_0.2.5.1                
[109] DirichletMultinomial_1.35.0   stringr_1.4.0                
[111] distributional_0.2.2          callr_3.7.0                  
[113] digest_0.6.27                 vegan_2.5-7                  
[115] graph_1.71.2                  rmarkdown_2.9                
[117] tidytree_0.3.4                DelayedMatrixStats_1.15.0    
[119] curl_4.3.2                    shiny_1.6.0                  
[121] lifecycle_1.0.0               nlme_3.1-152                 
[123] jsonlite_1.7.2                BiocNeighbors_1.11.0         
[125] desc_1.3.0                    CodeDepends_0.6.5            
[127] viridisLite_0.4.0             fansi_0.5.0                  
[129] pillar_1.6.2                  lattice_0.20-44              
[131] KEGGREST_1.33.0               fastmap_1.1.0                
[133] httr_1.4.2                    pkgbuild_1.2.0               
[135] interactiveDisplayBase_1.31.2 glue_1.4.2                   
[137] remotes_2.4.0                 png_0.1-7                    
[139] BiocVersion_3.14.0            bit_4.0.4                    
[141] stringi_1.7.3                 sass_0.4.0                   
[143] blob_1.2.2                    BiocSingular_1.9.1           
[145] AnnotationHub_3.1.4           memoise_2.0.0                
[147] dplyr_1.0.7                   irlba_2.3.3                  
[149] ape_5.5                      
```
</div>
