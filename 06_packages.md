# Packages {#packages}

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

<img src="general/figures/mia_logo.png" width="100" alt="mia logo" align="right" style="margin: 0 1em 0 1em" />

## Package installation

Several R packages provide methods for the analysis and manipulation
of `SummarizedExperiment` and related data containers. One of these is
`mia`. The installation for this and other packages has the following
procedure.

Stable Bioconductor release version can be installed with:


```r
BiocManager::install("microbiome/mia")
```

Biocondcuctor development version requires the installation of the
latest R beta version, and this is primarily recommended for those who
already have solid experience with R/Bioconductor and need access to
the latest experimental updates.


```r
BiocManager::install("microbiome/mia", version="devel")
```

The bleeding edge (and potentially unstable) development version lives
in Github:


```r
devtools::install_github("microbiome/mia")
```



## Some available packages

Some of the R packages supporting the _TreeSummarizedExperiment_ framework include:

### mia family of methods

- [mia](microbiome.github.io/mia) : generic microbiome analysis tools   
- [miaViz](microbiome.github.io/miaViz) : microbiome data visualization
- [miaSim](microbiome.github.io/miaSim) : microbiome data simulation
- [miaTime](microbiome.github.io/miaTime) : microbiome time series analysis

### Tree-based methods {sub-tree-methods}

- [philr](http://bioconductor.org/packages/devel/bioc/html/philr.html) (external, @Silverman2017)
- [mia](microbiome.github.io/mia): Microbiome analysis tools [@R-mia]
- [miaViz](microbiome.github.io/miaViz): Microbiome analysis specific visualization [@Ernst2022]
- [miaSim](microbiome.github.io/miaSim): Microbiome data simulations [@Simsek2021]
- [miaTime](microbiome.github.io/miaTime): Microbiome time series analysis [@Lahti2021]

### Differential abundance {sub-diff-abund}

- [benchdamic](https://bioconductor.org/packages/release/bioc/vignettes/benchdamic/inst/doc/intro.html) for benchmarking differential abundance methods
- [ANCOMBC](https://bioconductor.org/packages/devel/bioc/html/ANCOMBC.html) for differential abundance analysis

### Manipulation {sub-manipulation}

- [MicrobiotaProcess](https://bioconductor.org/packages/release/bioc/html/MicrobiotaProcess.html) for analyzing microbiome and other ecological data within the tidy framework

### Data

- [curatedMetagenomicData](https://bioconductor.org/packages/release/data/experiment/html/curatedMetagenomicData.html) a large collection of curated human microbiome data sets
- [microbiomeDataSets](https://bioconductor.org/packages/release/data/experiment/html/microbiomeDataSets.html) microbiome demo data sets

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
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] BiocStyle_2.24.0 rebook_1.6.0    

loaded via a namespace (and not attached):
 [1] bookdown_0.32       dir.expiry_1.4.0    codetools_0.2-19   
 [4] XML_3.99-0.13       digest_0.6.31       stats4_4.2.1       
 [7] evaluate_0.20       graph_1.74.0        rlang_1.0.6        
[10] cli_3.6.0           filelock_1.0.2      rmarkdown_2.20     
[13] tools_4.2.1         xfun_0.37           yaml_2.3.7         
[16] fastmap_1.1.1       compiler_4.2.1      BiocGenerics_0.44.0
[19] BiocManager_1.30.20 htmltools_0.5.4     CodeDepends_0.6.5  
[22] knitr_1.42         
```
</div>
