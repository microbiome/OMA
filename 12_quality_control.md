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

This chapter focuses on the exploration and quality control of
microbiome data and establishes commonly used descriptors of a
microbiome. The main difference to quality control is that the
exploration assumes that technical aspects of the dataset have been
investigated to your satisfaction. Generally speaking, at this point
you should be quite certain that the dataset does not suffer from
severe technical biases, or you should at least be aware of potential
problems.

In reality you might need to go back and forth between QC and exploration, 
since you discover through exploration of your dataset technical aspects you 
need to check.


```r
library(mia)
```


## Abundance

One initial approach for exploring data is by visualizing abundance. `miaViz` offers
the function `plotAbundanceDensity` where most abundant taxa can be plotted 
including several options.

In the following few demonstrations are shown, using the 
[@Lahti2014] dataset.

A Jitter plot based on relative abundance data, similar to the one presented at 
[@Salosensaari2021] supplementary figure 1, could be visualized as follows: 


```r
# Loading data
library(microbiomeDataSets)
tse <- atlas1006()

# Counts relative abundances
tse <- transformSamples(tse, method = "relabundance")

library(miaViz)
plotAbundanceDensity(tse, layout = "jitter", abund_values = "relabundance",
                     n = 40, point_size=1, point_shape=19, point_alpha=0.1) + 
    scale_x_log10(label=scales::percent)
```

<img src="12_quality_control_files/figure-html/unnamed-chunk-2-1.png" width="672" />

For instance, relative abundance values for the top 5 taxa can be
visualized as a density plot over a log scaled axis, using
"nationality" as an overlaying information:


```r
plotAbundanceDensity(tse, layout = "density", abund_values = "relabundance",
                     n = 5, colour_by="nationality", point_alpha=1/10 ) +
    scale_x_log10()
```

<img src="12_quality_control_files/figure-html/unnamed-chunk-3-1.png" width="672" />



## Prevalence

Prevalence is a measurement which describes in how many samples certain
microbes were detected.

Investigating the prevalence of microbes allows you either to focus on changes,
which pertain to most of the samples, or to focus on less often found microbes,
which are nonetheless abundantly found in a small number of samples.

On the raw data, the population prevalence (frequency) at a 1% relative
abundance threshold (`detection = 1/100` and `as_relative = TRUE`), can look
like this. The low prevalence in this example can be explained by rather
different sample types as well as the in-depth nature of the data.


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

The `detection` and `as_relative` can also be used to access, how many samples
do pass a threshold for raw counts. Here the population prevalence (frequency) 
at the absolute abundance threshold (`as_relative = FALSE`) at read count 1
(`detection = 1`) is accessed.


```r
head(getPrevalence(tse, detection = 1, sort = TRUE, abund_values = "counts",
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

Note that, if the output should be used for subsetting or storing the data in 
the `rowData`, set `sort = FALSE`.

### Prevalent microbiota analysis

To investigate the microbiome data at a selected taxonomic level, two 
approaches are available.

First the data can be agglomerated to the taxonomic level and `getPrevalence` 
be used on the result.


```r
altExp(tse,"Phylum") <- agglomerateByRank(tse, "Phylum")
head(getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = TRUE,
                   abund_values = "counts", as_relative = TRUE))
```

```
##      Firmicutes   Bacteroidetes  Actinobacteria  Proteobacteria Verrucomicrobia 
##       1.0000000       0.9852302       0.4821894       0.2988705       0.1277150 
##   Cyanobacteria 
##       0.0008688
```

Alternatively, the `rank` argument can be set to perform the agglomeration on
the fly.


```r
altExp(tse,"Phylum") <- agglomerateByRank(tse, "Phylum")
head(getPrevalence(tse, rank = "Phylum", detection = 1/100, sort = TRUE,
                   abund_values = "counts", as_relative = TRUE))
```

```
##      Firmicutes   Bacteroidetes  Actinobacteria  Proteobacteria Verrucomicrobia 
##       1.0000000       0.9852302       0.4821894       0.2988705       0.1277150 
##   Cyanobacteria 
##       0.0008688
```

The difference in the naming scheme is that, by default, `na.rm = TRUE` is used
for agglomeration in `getPrevalence`, whereas the default for 
`agglomerateByRank` is `FALSE` to prevent accidental data loss.

If you only need the names of the prevalent taxa, `getPrevalentTaxa` is
available. This returns the taxa that exceed the given prevalence and detection
thresholds.


```r
getPrevalentTaxa(tse, detection = 0, prevalence = 50/100)
prev <- getPrevalentTaxa(tse, detection = 0, prevalence = 50/100,
                         rank = "Phylum", sort = TRUE)
prev
```

Note that the `detection` and `prevalence` thresholds are not the same, since
`detection` can be applied to relative counts or absolute counts depending on 
whether `as_relative` is set `TRUE` or `FALSE`

TODO
See also related functions for the analysis of rare and variable taxa
(rareMembers; rareAbundance; lowAbundance). 

### Plotting prevalence

To plot the prevalence, the data is first added to the `rowData`.


```r
rowData(altExp(tse,"Phylum"))$prevalence <- 
    getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = FALSE,
                  abund_values = "counts", as_relative = TRUE)
```

Then it can be plotted via the plotting functions from the `scater` package.
 

```r
library(scater)
plotRowData(altExp(tse,"Phylum"), "prevalence", colour_by = "Phylum")
```

<img src="12_quality_control_files/figure-html/unnamed-chunk-7-1.png" width="672" />

Additionally, the prevalence can be plotted on the taxonomic tree using the
`miaViz` package.


```r
altExps(tse) <- splitByRanks(tse)
altExps(tse) <-
   lapply(altExps(tse),
          function(y){
              rowData(y)$prevalence <- 
                  getPrevalence(y, detection = 1/100, sort = FALSE,
                                abund_values = "counts", as_relative = TRUE)
              y
          })
top_phyla <- getTopTaxa(altExp(tse,"Phylum"),
                        method="prevalence",
                        top=5L,
                        abund_values="counts")
top_phyla_mean <- getTopTaxa(altExp(tse,"Phylum"),
                             method="mean",
                             top=5L,
                             abund_values="counts")
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse)[1:6])
x <- addTaxonomyTree(x)
```
 
After some preparation the data is assembled and can be plotted via 
`plotRowTree`.


```r
library(miaViz)
plotRowTree(x[rowData(x)$Phylum %in% top_phyla,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
```

<div class="figure">
<img src="12_quality_control_files/figure-html/plot-prev-prev-1.png" alt="Prevalence of top phyla as judged by prevalence" width="672" />
<p class="caption">(\#fig:plot-prev-prev)Prevalence of top phyla as judged by prevalence</p>
</div>

```r
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
```

<div class="figure">
<img src="12_quality_control_files/figure-html/plot-prev-mean-1.png" alt="Prevalence of top phyla as judged by mean abundance" width="672" />
<p class="caption">(\#fig:plot-prev-mean)Prevalence of top phyla as judged by mean abundance</p>
</div>





## Quality control


```r
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns 
```


### Top taxa  

The `getTopTaxa` can be used for identifying top taxa in the data.   

```r
top_features <- getTopTaxa(tse, method="median", top=10)
tax_data <- rowData(tse)[top_features,taxonomyRanks(tse)]
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


### Library size   

The total counts/sample can be calculated using the
`perCellQCMetrics`/`addPerCellQC` from the `scater` package. The former one 
just calculates the values, whereas the latter one directly adds them to the
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

<div class="figure">
<img src="12_quality_control_files/figure-html/plot-viz-lib-size-1-1.png" alt="Library size distribution." width="768" />
<p class="caption">(\#fig:plot-viz-lib-size-1)Library size distribution.</p>
</div>

Library sizes - and other variables from `colData` - can be also visualized by using 
specified function called `plotColData`.


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

<div class="figure">
<img src="12_quality_control_files/figure-html/plot-viz-lib-size-2-1.png" alt="Library sizes per sample." width="768" />
<p class="caption">(\#fig:plot-viz-lib-size-2)Library sizes per sample.</p>
</div>


```r
plotColData(tse,"sum","SampleType", colour_by = "SampleType") + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))
```

<div class="figure">
<img src="12_quality_control_files/figure-html/plot-viz-lib-size-3-1.png" alt="Library sizes per sample type." width="768" />
<p class="caption">(\#fig:plot-viz-lib-size-3)Library sizes per sample type.</p>
</div>


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
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] patchwork_1.1.1                dplyr_1.0.7                   
 [3] scater_1.22.0                  scuttle_1.4.0                 
 [5] miaViz_1.3.2                   ggraph_2.0.5                  
 [7] ggplot2_3.3.5                  microbiomeDataSets_1.1.5      
 [9] mia_1.3.13                     MultiAssayExperiment_1.20.0   
[11] TreeSummarizedExperiment_2.1.4 Biostrings_2.62.0             
[13] XVector_0.34.0                 SingleCellExperiment_1.16.0   
[15] SummarizedExperiment_1.24.0    Biobase_2.54.0                
[17] GenomicRanges_1.46.1           GenomeInfoDb_1.30.0           
[19] IRanges_2.28.0                 S4Vectors_0.32.3              
[21] BiocGenerics_0.40.0            MatrixGenerics_1.6.0          
[23] matrixStats_0.61.0-9001        BiocStyle_2.22.0              
[25] rebook_1.4.0                  

loaded via a namespace (and not attached):
  [1] AnnotationHub_3.2.0           BiocFileCache_2.2.0          
  [3] igraph_1.2.11                 plyr_1.8.6                   
  [5] lazyeval_0.2.2                splines_4.1.2                
  [7] BiocParallel_1.28.3           digest_0.6.29                
  [9] yulab.utils_0.0.4             htmltools_0.5.2              
 [11] viridis_0.6.2                 fansi_0.5.0                  
 [13] magrittr_2.0.1                memoise_2.0.1                
 [15] ScaledMatrix_1.2.0            cluster_2.1.2                
 [17] DECIPHER_2.22.0               graphlayouts_0.8.0           
 [19] colorspace_2.0-2              blob_1.2.2                   
 [21] rappdirs_0.3.3                ggrepel_0.9.1                
 [23] xfun_0.29                     crayon_1.4.2                 
 [25] RCurl_1.98-1.5                jsonlite_1.7.2               
 [27] graph_1.72.0                  ape_5.6                      
 [29] glue_1.6.0                    polyclip_1.10-0              
 [31] gtable_0.3.0                  zlibbioc_1.40.0              
 [33] DelayedArray_0.20.0           BiocSingular_1.10.0          
 [35] scales_1.1.1                  DBI_1.1.2                    
 [37] Rcpp_1.0.7                    viridisLite_0.4.0            
 [39] xtable_1.8-4                  decontam_1.14.0              
 [41] gridGraphics_0.5-1            tidytree_0.3.6               
 [43] bit_4.0.4                     rsvd_1.0.5                   
 [45] httr_1.4.2                    dir.expiry_1.2.0             
 [47] ellipsis_0.3.2                farver_2.1.0                 
 [49] pkgconfig_2.0.3               XML_3.99-0.8                 
 [51] CodeDepends_0.6.5             sass_0.4.0                   
 [53] dbplyr_2.1.1                  utf8_1.2.2                   
 [55] labeling_0.4.2                ggplotify_0.1.0              
 [57] tidyselect_1.1.1              rlang_0.4.12                 
 [59] reshape2_1.4.4                later_1.3.0                  
 [61] AnnotationDbi_1.56.2          munsell_0.5.0                
 [63] BiocVersion_3.14.0            tools_4.1.2                  
 [65] cachem_1.0.6                  DirichletMultinomial_1.36.0  
 [67] generics_0.1.1                RSQLite_2.2.9                
 [69] ExperimentHub_2.2.0           evaluate_0.14                
 [71] stringr_1.4.0                 fastmap_1.1.0                
 [73] yaml_2.2.1                    ggtree_3.2.1                 
 [75] knitr_1.37                    bit64_4.0.5                  
 [77] tidygraph_1.2.0               purrr_0.3.4                  
 [79] KEGGREST_1.34.0               nlme_3.1-153                 
 [81] sparseMatrixStats_1.6.0       mime_0.12                    
 [83] aplot_0.1.1                   compiler_4.1.2               
 [85] beeswarm_0.4.0                filelock_1.0.2               
 [87] curl_4.3.2                    png_0.1-7                    
 [89] interactiveDisplayBase_1.32.0 treeio_1.18.1                
 [91] tweenr_1.0.2                  tibble_3.1.6                 
 [93] bslib_0.3.1                   stringi_1.7.6                
 [95] highr_0.9                     lattice_0.20-45              
 [97] Matrix_1.4-0                  vegan_2.5-7                  
 [99] permute_0.9-5                 vctrs_0.3.8                  
[101] pillar_1.6.4                  lifecycle_1.0.1              
[103] BiocManager_1.30.16           jquerylib_0.1.4              
[105] BiocNeighbors_1.12.0          cowplot_1.1.1                
[107] bitops_1.0-7                  irlba_2.3.5                  
[109] httpuv_1.6.5                  R6_2.5.1                     
[111] bookdown_0.24                 promises_1.2.0.1             
[113] gridExtra_2.3                 vipor_0.4.5                  
[115] codetools_0.2-18              MASS_7.3-54                  
[117] assertthat_0.2.1              withr_2.4.3                  
[119] GenomeInfoDbData_1.2.7        mgcv_1.8-38                  
[121] parallel_4.1.2                ggfun_0.0.4                  
[123] grid_4.1.2                    beachmat_2.10.0              
[125] tidyr_1.1.4                   rmarkdown_2.11               
[127] DelayedMatrixStats_1.16.0     ggnewscale_0.4.5             
[129] ggforce_0.3.3                 shiny_1.7.1                  
[131] ggbeeswarm_0.6.0             
```
</div>

