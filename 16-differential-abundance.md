# Differential abundance {#differential-abundance}

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


## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R Under development (unstable) (2021-04-08 r80148)
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
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] BiocStyle_2.19.2 rebook_1.1.20   

loaded via a namespace (and not attached):
 [1] graph_1.69.0        knitr_1.31          magrittr_2.0.1     
 [4] BiocGenerics_0.37.1 R6_2.5.0            rlang_0.4.10       
 [7] stringr_1.4.0       tools_4.1.0         parallel_4.1.0     
[10] xfun_0.22           jquerylib_0.1.3     htmltools_0.5.1.1  
[13] CodeDepends_0.6.5   yaml_2.2.1          digest_0.6.27      
[16] bookdown_0.21       dir.expiry_0.99.4   BiocManager_1.30.12
[19] codetools_0.2-18    sass_0.3.1          evaluate_0.14      
[22] rmarkdown_2.7       stringi_1.5.3       compiler_4.1.0     
[25] bslib_0.2.4         filelock_1.0.2      XML_3.99-0.6       
[28] stats4_4.1.0        jsonlite_1.7.2     
```
</div>
