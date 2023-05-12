# Exploration and quality Control {#quality-control}


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

This chapter focuses on the quality control and exploration of
microbiome data and establishes commonly used descriptive
summaries. Familiarizing with the peculiarities of a given data set is
the essential basis for any data analysis and model building.

The dataset should not suffer from severe technical biases, and you
should at least be aware of potential challenges, such as outliers,
biases, unexpected patterns and so forth. Standard summaries and
visualizations can help, and the rest comes with experience. The
exploration and quality control can be iterative processes.



```r
library(mia)
```


## Abundance

Abundance visualization is an important data exploration
approach. `miaViz` offers the function `plotAbundanceDensity` to plot
the most abundant taxa with several options.

Next, a few demonstrations are shown, using the [@Lahti2014]
dataset. A Jitter plot based on relative abundance data, similar to
the one presented at [@Salosensaari2021] supplementary figure 1, could
be visualized as follows:


```r
# Load example data
library(miaTime)
library(miaViz)
data(hitchip1006)
tse <- hitchip1006

# Add relative abundances
tse <- transformCounts(tse, MARGIN = "samples", method = "relabundance")

# Use argument names
# assay.type / assay.type / assay.type
# depending on the mia package version
plotAbundanceDensity(tse, layout = "jitter", assay.type = "relabundance",
                     n = 40, point_size=1, point_shape=19, point_alpha=0.1) + 
                     scale_x_log10(label=scales::percent)
```

![](12_quality_control_files/figure-latex/unnamed-chunk-2-1.pdf)<!-- --> 

The relative abundance values for the top-5 taxonomic features can be
visualized as a density plot over a log scaled axis, with
"nationality" indicated by colors:


```r
plotAbundanceDensity(tse, layout = "density", assay.type = "relabundance",
                     n = 5, colour_by="nationality", point_alpha=1/10) +
    scale_x_log10()
```

![](12_quality_control_files/figure-latex/unnamed-chunk-3-1.pdf)<!-- --> 



## Prevalence

Prevalence quantifies the frequency of samples where certain microbes
were detected (above a given detection threshold). The prevalence can
be given as sample size (N) or percentage (unit interval).

Investigating prevalence allows you either to focus on changes which
pertain to the majority of the samples, or identify rare microbes,
which may be _conditionally abundant_ in a small number of samples.

The population prevalence (frequency) at a 1% relative abundance
threshold (`detection = 1/100` and `as_relative = TRUE`), can look
like this. 


```r
head(getPrevalence(tse, detection = 1/100, sort = TRUE, as_relative = TRUE))
```

```
## Faecalibacterium prausnitzii et rel.           Ruminococcus obeum et rel. 
##                               0.9522                               0.9140 
##   Oscillospira guillermondii et rel.        Clostridium symbiosum et rel. 
##                               0.8801                               0.8714 
##     Subdoligranulum variable at rel.     Clostridium orbiscindens et rel. 
##                               0.8358                               0.8315
```

The function arguments `detection` and `as_relative` can also be used
to access, how many samples do pass a threshold for raw counts. Here,
the population prevalence (frequency) at the absolute abundance
threshold (`as_relative = FALSE`) at read count 1 (`detection = 1`) is
accessed.


```r
head(getPrevalence(tse, detection = 1, sort = TRUE, assay.type = "counts",
                   as_relative = FALSE))
```

```
##            Uncultured Mollicutes      Uncultured Clostridiales II 
##                                1                                1 
##       Uncultured Clostridiales I               Tannerella et rel. 
##                                1                                1 
##   Sutterella wadsworthia et rel. Subdoligranulum variable at rel. 
##                                1                                1
```

If the output should be used for subsetting or storing the data in the
`rowData`, set `sort = FALSE`.


### Prevalence analysis

To investigate microbiome prevalence at a selected taxonomic level, two 
approaches are available.

First the data can be agglomerated to the taxonomic level and `getPrevalence` 
applied on the resulting object.


```r
# Agglomerate taxa abundances to Phylum level, and add the new table
# to the altExp slot
altExp(tse,"Phylum") <- agglomerateByRank(tse, "Phylum")
# Check prevalence for the Phylum abundance table from the altExp slot
head(getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = TRUE,
                   assay.type = "counts", as_relative = TRUE))
```

```
##      Firmicutes   Bacteroidetes  Actinobacteria  Proteobacteria Verrucomicrobia 
##       1.0000000       0.9852302       0.4821894       0.2988705       0.1277150 
##   Cyanobacteria 
##       0.0008688
```


Alternatively, the `rank` argument could be set to perform the
agglomeration on the fly.


```r
head(getPrevalence(tse, rank = "Phylum", detection = 1/100, sort = TRUE,
                   assay.type = "counts", as_relative = TRUE))
```

```
##      Firmicutes   Bacteroidetes  Actinobacteria  Proteobacteria Verrucomicrobia 
##       1.0000000       0.9852302       0.4821894       0.2988705       0.1277150 
##   Cyanobacteria 
##       0.0008688
```

Note that, by default, `na.rm = TRUE` is used for agglomeration in
`getPrevalence`, whereas the default for `agglomerateByRank` is
`FALSE` to prevent accidental data loss.

If you only need the names of the prevalent taxa, `getPrevalentTaxa`
is available. This returns the taxa that exceed the given prevalence
and detection thresholds.


```r
getPrevalentTaxa(tse, detection = 0, prevalence = 50/100)
prev <- getPrevalentTaxa(tse, detection = 0, prevalence = 50/100,
                         rank = "Phylum", sort = TRUE)
prev
```

Note that the `detection` and `prevalence` thresholds are not the same, since
`detection` can be applied to relative counts or absolute counts depending on 
whether `as_relative` is set `TRUE` or `FALSE`


The function ‘getPrevalentAbundance’ can be used to check the total
relative abundance of the prevalent taxa (between 0 and 1).


### Rare taxa

Related functions are available for the analysis of rare taxa
(`rareMembers`; `rareAbundance`; `lowAbundance`, `getRareTaxa`,
`subsetByRareTaxa`).


### Plotting prevalence

To plot the prevalence, add the prevalence of each taxon to
`rowData`. Here, we are analysing the Phylum level abundances, which
are stored in the `altExp` slot.


```r
rowData(altExp(tse,"Phylum"))$prevalence <- 
    getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = FALSE,
                  assay.type = "counts", as_relative = TRUE)
```

The prevalences can then be plotted using the plotting functions from
the `scater` package.
 

```r
library(scater)
plotRowData(altExp(tse,"Phylum"), "prevalence", colour_by = "Phylum")
```

![](12_quality_control_files/figure-latex/unnamed-chunk-7-1.pdf)<!-- --> 

The prevalence can also be visualized on the taxonomic tree with the
`miaViz` package.


```r
altExps(tse) <- splitByRanks(tse)
altExps(tse) <-
   lapply(altExps(tse),
          function(y){
              rowData(y)$prevalence <- 
                  getPrevalence(y, detection = 1/100, sort = FALSE,
                                assay.type = "counts", as_relative = TRUE)
              y
          })
top_phyla <- getTopTaxa(altExp(tse,"Phylum"),
                        method="prevalence",
                        top=5L,
                        assay.type="counts")
top_phyla_mean <- getTopTaxa(altExp(tse,"Phylum"),
                             method="mean",
                             top=5L,
                             assay.type="counts")
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse)[1:6])
x <- addTaxonomyTree(x)
```
 
After some preparation, the data is assembled and can be plotted with
`plotRowTree`.


```r
library(miaViz)
plotRowTree(x[rowData(x)$Phylum %in% top_phyla,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
```

![(\#fig:plot-prev-prev)Prevalence of top phyla as judged by prevalence](12_quality_control_files/figure-latex/plot-prev-prev-1.pdf) 



```r
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
```

![(\#fig:plot-prev-mean)Prevalence of top phyla as judged by mean abundance](12_quality_control_files/figure-latex/plot-prev-mean-1.pdf) 

## Quality control {#qc}

Next, let us load the `GlobalPatterns` data set to illustrate standard
microbiome data summaries.


```r
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns 
```


### Top taxa  

The `getTopTaxa` identifies top taxa in the data.   


```r
# Pick the top taxa
top_features <- getTopTaxa(tse, method="median", top=10)

# Check the information for these
rowData(tse)[top_features, taxonomyRanks(tse)]
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


### Library size / read count  

The total counts/sample can be calculated using `perCellQCMetrics`/`addPerCellQC` from the `scater` package. The former one
just calculates the values, whereas the latter one directly adds them to
`colData`.


```r
library(scater)
perCellQCMetrics(tse)
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
tse <- addPerCellQC(tse)
colData(tse)
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

The distribution of calculated library sizes can be visualized as a
histogram (left), or by sorting the samples by library size (right).


```r
library(ggplot2)

p1 <- ggplot(as.data.frame(colData(tse))) +
        geom_histogram(aes(x = sum), color = "black", fill = "gray", bins = 30) +
        labs(x = "Library size", y = "Frequency (n)") + 
        # scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
        # labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), # Removes the grid
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) # Adds y-axis

library(dplyr)
df <- as.data.frame(colData(tse)) %>%
        arrange(sum) %>%
        mutate(index = 1:n())
p2 <- ggplot(df, aes(y = index, x = sum/1e6)) +
        geom_point() +	
        labs(x = "Library size (million reads)", y = "Sample index") +	
        theme_bw() +
        theme(panel.grid.major = element_blank(), # Removes the grid
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) # Adds y-axis

library(patchwork)
p1 + p2
```

![(\#fig:plot-viz-lib-size-1)Library size distribution.](12_quality_control_files/figure-latex/plot-viz-lib-size-1-1.pdf) 

Library sizes other variables from `colData` can be
visualized by using specified function called `plotColData`.


```r
library(ggplot2)
# Sort samples by read count, order the factor levels, and store back to tse as DataFrame
# TODO: plotColData could include an option for sorting samples based on colData variables
colData(tse) <- as.data.frame(colData(tse)) %>%
                 arrange(X.SampleID) %>%
        	 mutate(X.SampleID = factor(X.SampleID, levels=X.SampleID)) %>%
		 DataFrame
plotColData(tse,"sum","X.SampleID", colour_by = "SampleType") + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    labs(y = "Library size (N)", x = "Sample ID") 	    
```

![(\#fig:plot-viz-lib-size-2)Library sizes per sample.](12_quality_control_files/figure-latex/plot-viz-lib-size-2-1.pdf) 


```r
plotColData(tse,"sum","SampleType", colour_by = "SampleType") + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))
```

![(\#fig:plot-viz-lib-size-3)Library sizes per sample type.](12_quality_control_files/figure-latex/plot-viz-lib-size-3-1.pdf) 

In addition, data can be rarefied with
[subsampleCounts](https://microbiome.github.io/mia/reference/subsampleCounts.html),
which normalises the samples to an equal number of reads. However,
this practice has been discouraged for the analysis of differentially
abundant microorganisms (see [@mcmurdie2014waste]).
  

### Contaminant sequences

Samples might be contaminated with exogenous sequences. The impact of
each contaminant can be estimated based on their frequencies and
concentrations across the samples.

The following [decontam
functions](https://microbiome.github.io/mia/reference/isContaminant.html)
are based on the [@davis2018simple] and support such functionality:

* `isContaminant`, `isNotContaminant`
* `addContaminantQC`, `addNotContaminantQC`


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
 [1] patchwork_1.1.2                dplyr_1.1.2                   
 [3] scater_1.28.0                  scuttle_1.10.1                
 [5] miaViz_1.9.1                   ggraph_2.1.0                  
 [7] ggplot2_3.4.2                  miaTime_0.1.21                
 [9] mia_1.9.2                      MultiAssayExperiment_1.26.0   
[11] TreeSummarizedExperiment_2.1.4 Biostrings_2.68.0             
[13] XVector_0.40.0                 SingleCellExperiment_1.22.0   
[15] SummarizedExperiment_1.30.1    Biobase_2.60.0                
[17] GenomicRanges_1.52.0           GenomeInfoDb_1.36.0           
[19] IRanges_2.34.0                 S4Vectors_0.38.1              
[21] BiocGenerics_0.46.0            MatrixGenerics_1.12.0         
[23] matrixStats_0.63.0-9003        BiocStyle_2.28.0              
[25] rebook_1.9.0                  

loaded via a namespace (and not attached):
  [1] jsonlite_1.8.4              CodeDepends_0.6.5          
  [3] magrittr_2.0.3              ggbeeswarm_0.7.2           
  [5] farver_2.1.1                rmarkdown_2.21             
  [7] zlibbioc_1.46.0             vctrs_0.6.2                
  [9] memoise_2.0.1               DelayedMatrixStats_1.22.0  
 [11] RCurl_1.98-1.12             ggtree_3.8.0               
 [13] htmltools_0.5.5             S4Arrays_1.0.1             
 [15] BiocNeighbors_1.18.0        gridGraphics_0.5-1         
 [17] plyr_1.8.8                  DECIPHER_2.28.0            
 [19] cachem_1.0.8                igraph_1.4.2               
 [21] lifecycle_1.0.3             pkgconfig_2.0.3            
 [23] rsvd_1.0.5                  Matrix_1.5-4               
 [25] R6_2.5.1                    fastmap_1.1.1              
 [27] GenomeInfoDbData_1.2.10     digest_0.6.31              
 [29] aplot_0.1.10                colorspace_2.1-0           
 [31] ggnewscale_0.4.8            irlba_2.3.5.1              
 [33] RSQLite_2.3.1               vegan_2.6-4                
 [35] beachmat_2.16.0             labeling_0.4.2             
 [37] filelock_1.0.2              fansi_1.0.4                
 [39] polyclip_1.10-4             mgcv_1.8-42                
 [41] compiler_4.3.0              bit64_4.0.5                
 [43] withr_2.5.0                 BiocParallel_1.34.1        
 [45] viridis_0.6.3               DBI_1.1.3                  
 [47] highr_0.10                  ggforce_0.4.1              
 [49] MASS_7.3-60                 DelayedArray_0.26.2        
 [51] permute_0.9-7               tools_4.3.0                
 [53] vipor_0.4.5                 beeswarm_0.4.0             
 [55] ape_5.7-1                   glue_1.6.2                 
 [57] nlme_3.1-162                grid_4.3.0                 
 [59] cluster_2.1.4               reshape2_1.4.4             
 [61] generics_0.1.3              gtable_0.3.3               
 [63] tidyr_1.3.0                 BiocSingular_1.16.0        
 [65] tidygraph_1.2.3             ScaledMatrix_1.8.1         
 [67] utf8_1.2.3                  ggrepel_0.9.3              
 [69] pillar_1.9.0                stringr_1.5.0              
 [71] yulab.utils_0.0.6           splines_4.3.0              
 [73] tweenr_2.0.2                treeio_1.24.0              
 [75] lattice_0.21-8              bit_4.0.5                  
 [77] tidyselect_1.2.0            DirichletMultinomial_1.42.0
 [79] knitr_1.42                  gridExtra_2.3              
 [81] bookdown_0.34               xfun_0.39                  
 [83] graphlayouts_1.0.0          stringi_1.7.12             
 [85] lazyeval_0.2.2              ggfun_0.0.9                
 [87] yaml_2.3.7                  evaluate_0.21              
 [89] codetools_0.2-19            tibble_3.2.1               
 [91] BiocManager_1.30.20         graph_1.78.0               
 [93] ggplotify_0.1.0             cli_3.6.1                  
 [95] munsell_0.5.0               Rcpp_1.0.10                
 [97] dir.expiry_1.8.0            XML_3.99-0.14              
 [99] parallel_4.3.0              blob_1.2.4                 
[101] sparseMatrixStats_1.12.0    bitops_1.0-7               
[103] decontam_1.20.0             viridisLite_0.4.2          
[105] tidytree_0.4.2              scales_1.2.1               
[107] purrr_1.0.1                 crayon_1.5.2               
[109] rlang_1.1.1                 cowplot_1.1.1              
```
</div>

