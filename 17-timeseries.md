# Microbiome time series {#microbiome-timeseries}

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

This chapter focuses on the exploration of microbiome time series.


## Stability

## Tipping elements


## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.0.3 (2020-10-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04 LTS

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
[1] BiocStyle_2.18.1    rebook_1.0.0        BiocManager_1.30.10

loaded via a namespace (and not attached):
 [1] bookdown_0.21       codetools_0.2-18    XML_3.99-0.5       
 [4] ps_1.5.0            digest_0.6.27       stats4_4.0.3       
 [7] magrittr_2.0.1      evaluate_0.14       graph_1.68.0       
[10] rlang_0.4.9         stringi_1.5.3       callr_3.5.1        
[13] rmarkdown_2.5       tools_4.0.3         stringr_1.4.0      
[16] processx_3.4.5      parallel_4.0.3      xfun_0.19          
[19] yaml_2.2.1          compiler_4.0.3      BiocGenerics_0.36.0
[22] htmltools_0.5.0     CodeDepends_0.6.5   knitr_1.30         
```
</div>
