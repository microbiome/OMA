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

You can install all packages that are required to run every example in this book via the following command:


```r
source("https://raw.githubusercontent.com/microbiome/OMA/master/install_packages.R")
```


### Installing specific packages {#packages}

You can install R packages of your choice with the following command
line procedure.

**Bioconductor release version** is the most stable and tested version
but may miss some of the latest methods and updates. It can be
installed with:


```r
BiocManager::install("microbiome/mia")
```

**Bioconductor development version** requires the installation of the
latest R beta version. This is primarily recommended for those who
already have experience with R/Bioconductor and need access to the
latest updates.


```r
BiocManager::install("microbiome/mia", version="devel")
```

**Github development version** provides access to the latest but
potentially unstable features. This is useful when you want access to
all available tools.


```r
devtools::install_github("microbiome/mia")
```


## Package ecosystem 

Methods for the analysis and manipulation of
`(Tree)SummarizedExperiment` and `MultiAssayExperiment` data
containers are available through a number of R packages. Some of these
are listed below. If you know more tips on such packages, data
sources, or other resources, kindly [let us
know](https://microbiome.github.io) through the issues, pull requests,
or online channels.


### mia family of methods

- [mia](microbiome.github.io/mia): Microbiome analysis tools [@R-mia]
- [miaViz](microbiome.github.io/miaViz): Microbiome analysis specific visualization [@Ernst2022]
- [miaSim](microbiome.github.io/miaSim): Microbiome data simulations [@Simsek2021]
- [miaTime](microbiome.github.io/miaTime): Microbiome time series analysis [@Lahti2021]


### Tree-based methods {#sub-tree-methods}

- [philr](http://bioconductor.org/packages/devel/bioc/html/philr.html) (@Silverman2017)


### Differential abundance {#sub-diff-abund}

- [ANCOMBC](https://bioconductor.org/packages/devel/bioc/html/ANCOMBC.html) for differential abundance analysis
- [benchdamic](https://bioconductor.org/packages/release/bioc/vignettes/benchdamic/inst/doc/intro.html) for benchmarking differential abundance methods

### Manipulation {#sub-manipulation}

- [MicrobiotaProcess](https://bioconductor.org/packages/release/bioc/html/MicrobiotaProcess.html) for analyzing microbiome and other ecological data within the tidy framework


### Further options

- [Tools for Microbiome
  Analysis](https://microsud.github.io/Tools-Microbiome-Analysis/)
  site listed over 130 R packages for microbiome data science in
  2023. Many of these do not directly support the data containers used
  in this book but can be used with minor conversions between formats.

## Data ecosystem

Section \@ref(example-data) links to various open microbiome data resources that support this framework.



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
 [1] bookdown_0.33       dir.expiry_1.4.0    codetools_0.2-19   
 [4] XML_3.99-0.14       digest_0.6.31       stats4_4.2.1       
 [7] evaluate_0.20       graph_1.74.0        rlang_1.1.0        
[10] cli_3.6.1           filelock_1.0.2      rmarkdown_2.21     
[13] tools_4.2.1         xfun_0.38           yaml_2.3.7         
[16] fastmap_1.1.1       compiler_4.2.1      BiocGenerics_0.44.0
[19] BiocManager_1.30.20 htmltools_0.5.5     CodeDepends_0.6.5  
[22] knitr_1.42         
```
</div>
