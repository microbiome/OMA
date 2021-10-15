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


## Package installation

Several R packages are now available providing methods for the
analysis and manipulation of `SummarizedExperiment` data and its
derivatives. One of these is `mia`.

Install the Biocondcuctor development version with:


```r
BiocManager::install("microbiome/mia", version="devel")
```


## Available packages

Some of the R packages supporting the framework include:

- [mia](microbiome.github.io/mia) : Microbiome analysis tools   
- [miaViz](microbiome.github.io/miaViz) : Microbiome analysis specific visualization
- [miaSim](microbiome.github.io/miaSim) : Microbiome data simulations
- [miaTime](microbiome.github.io/miaTime) : Microbiome time series analysis
- [philr](http://bioconductor.org/packages/devel/bioc/html/philr.html) (external)



## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.1.1 (2021-08-10)
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
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] BiocStyle_2.21.4 rebook_1.3.1    

loaded via a namespace (and not attached):
 [1] graph_1.71.2        knitr_1.36          magrittr_2.0.1     
 [4] BiocGenerics_0.39.2 R6_2.5.1            rlang_0.4.11       
 [7] fastmap_1.1.0       stringr_1.4.0       tools_4.1.1        
[10] xfun_0.26           jquerylib_0.1.4     htmltools_0.5.2    
[13] CodeDepends_0.6.5   yaml_2.2.1          digest_0.6.28      
[16] bookdown_0.24       dir.expiry_1.1.0    BiocManager_1.30.16
[19] codetools_0.2-18    sass_0.4.0          evaluate_0.14      
[22] rmarkdown_2.11      stringi_1.7.5       compiler_4.1.1     
[25] bslib_0.3.1         filelock_1.0.2      stats4_4.1.1       
[28] XML_3.99-0.8        jsonlite_1.7.2     
```
</div>
