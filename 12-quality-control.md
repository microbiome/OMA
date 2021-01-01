# Quality Control {#quality-control}

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
data("GlobalPatterns")
se <- GlobalPatterns 
```

## Get top taxa and taxonomy   

### Features  

The `getTopTaxa` can be used for identifying top taxa in the data.   

```r
top_features <- getTopTaxa(se, method="median", top=10)
tax_data <- rowData(se)[top_features,taxonomyRanks(se)]
tax_data
```

```
## DataFrame with 10 rows and 7 columns
##            Kingdom         Phylum               Class             Order
##        <character>    <character>         <character>       <character>
## 549656    Bacteria  Cyanobacteria         Chloroplast     Stramenopiles
## 331820    Bacteria  Bacteroidetes         Bacteroidia     Bacteroidales
## 317182    Bacteria  Cyanobacteria         Chloroplast     Stramenopiles
## 94166     Bacteria Proteobacteria Gammaproteobacteria    Pasteurellales
## 279599    Bacteria  Cyanobacteria    Nostocophycideae        Nostocales
## 158660    Bacteria  Bacteroidetes         Bacteroidia     Bacteroidales
## 329744    Bacteria Actinobacteria      Actinobacteria   Actinomycetales
## 326977    Bacteria Actinobacteria      Actinobacteria Bifidobacteriales
## 248140    Bacteria  Bacteroidetes         Bacteroidia     Bacteroidales
## 550960    Bacteria Proteobacteria Gammaproteobacteria Enterobacteriales
##                    Family           Genus                Species
##               <character>     <character>            <character>
## 549656                 NA              NA                     NA
## 331820     Bacteroidaceae     Bacteroides                     NA
## 317182                 NA              NA                     NA
## 94166     Pasteurellaceae     Haemophilus Haemophilusparainflu..
## 279599        Nostocaceae  Dolichospermum                     NA
## 158660     Bacteroidaceae     Bacteroides                     NA
## 329744             ACK-M1              NA                     NA
## 326977 Bifidobacteriaceae Bifidobacterium Bifidobacteriumadole..
## 248140     Bacteroidaceae     Bacteroides      Bacteroidescaccae
## 550960 Enterobacteriaceae     Providencia                     NA
```

## Library size   

The total counts/sample can be calculated using the
`perCellQCMetrics`/`addPerCellQC` from the `scater` package. The first one 
just calculates the values, whereas the latter one directly adds them to the
`colData`.


```r
library(scater)
perCellQCMetrics(se)
```

```
## DataFrame with 26 rows and 3 columns
##               sum  detected     total
##         <numeric> <numeric> <numeric>
## CL3        864077      6964    864077
## CC1       1135457      7679   1135457
## SV1        697509      5729    697509
## M31Fcsw   1543451      2667   1543451
## M11Fcsw   2076476      2574   2076476
## ...           ...       ...       ...
## TS28       937466      2679    937466
## TS29      1211071      2629   1211071
## Even1     1216137      4213   1216137
## Even2      971073      3130    971073
## Even3     1078241      2776   1078241
```

```r
se <- addPerCellQC(se)
colData(se)
```

```
## DataFrame with 26 rows and 10 columns
##         X.SampleID   Primer Final_Barcode Barcode_truncated_plus_T
##           <factor> <factor>      <factor>                 <factor>
## CL3        CL3      ILBC_01        AACGCA                   TGCGTT
## CC1        CC1      ILBC_02        AACTCG                   CGAGTT
## SV1        SV1      ILBC_03        AACTGT                   ACAGTT
## M31Fcsw    M31Fcsw  ILBC_04        AAGAGA                   TCTCTT
## M11Fcsw    M11Fcsw  ILBC_05        AAGCTG                   CAGCTT
## ...            ...      ...           ...                      ...
## TS28         TS28   ILBC_25        ACCAGA                   TCTGGT
## TS29         TS29   ILBC_26        ACCAGC                   GCTGGT
## Even1        Even1  ILBC_27        ACCGCA                   TGCGGT
## Even2        Even2  ILBC_28        ACCTCG                   CGAGGT
## Even3        Even3  ILBC_29        ACCTGT                   ACAGGT
##         Barcode_full_length SampleType
##                    <factor>   <factor>
## CL3             CTAGCGTGCGT      Soil 
## CC1             CATCGACGAGT      Soil 
## SV1             GTACGCACAGT      Soil 
## M31Fcsw         TCGACATCTCT      Feces
## M11Fcsw         CGACTGCAGCT      Feces
## ...                     ...        ...
## TS28            GCATCGTCTGG      Feces
## TS29            CTAGTCGCTGG      Feces
## Even1           TGACTCTGCGG      Mock 
## Even2           TCTGATCGAGG      Mock 
## Even3           AGAGAGACAGG      Mock 
##                                        Description       sum  detected
##                                           <factor> <numeric> <numeric>
## CL3     Calhoun South Carolina Pine soil, pH 4.9      864077      6964
## CC1     Cedar Creek Minnesota, grassland, pH 6.1     1135457      7679
## SV1     Sevilleta new Mexico, desert scrub, pH 8.3    697509      5729
## M31Fcsw M3, Day 1, fecal swab, whole body study      1543451      2667
## M11Fcsw M1, Day 1, fecal swab, whole body study      2076476      2574
## ...                                            ...       ...       ...
## TS28                                       Twin #1    937466      2679
## TS29                                       Twin #2   1211071      2629
## Even1                                      Even1     1216137      4213
## Even2                                      Even2      971073      3130
## Even3                                      Even3     1078241      2776
##             total
##         <numeric>
## CL3        864077
## CC1       1135457
## SV1        697509
## M31Fcsw   1543451
## M11Fcsw   2076476
## ...           ...
## TS28       937466
## TS29      1211071
## Even1     1216137
## Even2      971073
## Even3     1078241
```

The calulated library sizes can be visualized using the `plotColData` function.


```r
library(ggplot2)
plotColData(se,"sum","X.SampleID", colour_by = "X.SampleID") + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))
```

<div class="figure">
<img src="12-quality-control_files/figure-html/plot-viz-lib-size-1-1.png" alt="Library sizes per sample." width="768" />
<p class="caption">(\#fig:plot-viz-lib-size-1)Library sizes per sample.</p>
</div>


```r
plotColData(se,"sum","SampleType", colour_by = "SampleType") + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))
```

<div class="figure">
<img src="12-quality-control_files/figure-html/plot-viz-lib-size-2-1.png" alt="Library sizes per sample type." width="768" />
<p class="caption">(\#fig:plot-viz-lib-size-2)Library sizes per sample type.</p>
</div>


## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.0.3 (2020-10-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.1 LTS

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
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] scater_1.18.3                    ggplot2_3.3.3                   
 [3] mia_0.98.21                      MicrobiomeExperiment_0.99.0.9014
 [5] Biostrings_2.58.0                XVector_0.30.0                  
 [7] TreeSummarizedExperiment_1.6.2   SingleCellExperiment_1.12.0     
 [9] SummarizedExperiment_1.20.0      Biobase_2.50.0                  
[11] GenomicRanges_1.42.0             GenomeInfoDb_1.26.2             
[13] IRanges_2.24.1                   S4Vectors_0.28.1                
[15] BiocGenerics_0.36.0              MatrixGenerics_1.2.0            
[17] matrixStats_0.57.0               BiocStyle_2.18.1                
[19] rebook_1.0.0                     BiocManager_1.30.10             

loaded via a namespace (and not attached):
 [1] nlme_3.1-151                bitops_1.0-6               
 [3] DirichletMultinomial_1.32.0 tools_4.0.3                
 [5] R6_2.5.0                    irlba_2.3.3                
 [7] vegan_2.5-7                 vipor_0.4.5                
 [9] mgcv_1.8-33                 colorspace_2.0-0           
[11] permute_0.9-5               withr_2.3.0                
[13] tidyselect_1.1.0            gridExtra_2.3              
[15] processx_3.4.5              compiler_4.0.3             
[17] graph_1.68.0                BiocNeighbors_1.8.2        
[19] DelayedArray_0.16.0         labeling_0.4.2             
[21] bookdown_0.21               scales_1.1.1               
[23] callr_3.5.1                 stringr_1.4.0              
[25] digest_0.6.27               rmarkdown_2.6              
[27] pkgconfig_2.0.3             htmltools_0.5.0            
[29] sparseMatrixStats_1.2.0     highr_0.8                  
[31] rlang_0.4.10                DelayedMatrixStats_1.12.1  
[33] farver_2.0.3                generics_0.1.0             
[35] BiocParallel_1.24.1         dplyr_1.0.2                
[37] RCurl_1.98-1.2              magrittr_2.0.1             
[39] BiocSingular_1.6.0          GenomeInfoDbData_1.2.4     
[41] scuttle_1.0.4               Matrix_1.3-0               
[43] Rcpp_1.0.5                  ggbeeswarm_0.6.0           
[45] munsell_0.5.0               ape_5.4-1                  
[47] viridis_0.5.1               lifecycle_0.2.0            
[49] stringi_1.5.3               yaml_2.2.1                 
[51] MASS_7.3-53                 zlibbioc_1.36.0            
[53] grid_4.0.3                  crayon_1.3.4               
[55] lattice_0.20-41             cowplot_1.1.1              
[57] beachmat_2.6.4              splines_4.0.3              
[59] CodeDepends_0.6.5           knitr_1.30                 
[61] ps_1.5.0                    pillar_1.4.7               
[63] codetools_0.2-18            XML_3.99-0.5               
[65] glue_1.4.2                  evaluate_0.14              
[67] vctrs_0.3.6                 gtable_0.3.0               
[69] purrr_0.3.4                 tidyr_1.1.2                
[71] xfun_0.19                   rsvd_1.0.3                 
[73] viridisLite_0.3.0           tibble_3.0.4               
[75] beeswarm_0.2.3              cluster_2.1.0              
[77] ellipsis_0.3.1             
```
</div>
