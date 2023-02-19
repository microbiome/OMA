# (PART) Appendix {-}



# Training {#training}

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

The tutorial can support course teaching in R/Bioconductor.


## Code of Conduct {#coc}

We support the Bioconductor Code of Conduct. In short, the Bioconductor community values an open approach to science that promotes 

    * sharing of ideas, code, software and expertise
    * collaboration
    * diversity and inclusivity
    * a kind and welcoming environment
    * community contributions

For further details, we refer the reader to the full [Code of Conduct](https://bioconductor.github.io/bioc_coc_multilingual/).


## Checklist {#checklist}

### Setting up Your own computer

Setting up the system on your own computer is not required for the
course but it can be useful for later use. The required software:

* [R (the latest official release)](https://www.r-project.org/) 

* [RStudio](https://www.rstudio.com/products/rstudio/download/);
  choose "Rstudio Desktop" to download the latest version. RStudio is
  optional. Check the [Rstudio home page](https://www.rstudio.com/)
  for more information.

* Install and load essential R packages (see Section \@ref(packages))

* After a successful installation you can start with the
  case study examples in Section \@ref(exercises).


### Support and resources

 * Check additional reading tips and try out online material listed in
   Section \@ref(material).

 * **You can run the workflows by simply copy-pasting the examples.**
  You can then test further examples from this tutorial, modifying and
  applying these techniques to your own data.

 * Online support on installation and other matters, join us at
   [Gitter](https://gitter.im/microbiome/miaverse?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


### Study material {#material}

We encourage to familiarize with the material and test examples in advance.

 * [Orchestrating Microbiome Analysis with R/Bioconductor (OMA)](https://microbiome.github.io/OMA/)

 * Short online videos on [R/Bioconductor microbiome data science tools](https://www.youtube.com/playlist?list=PLjiXAZO27elAJEptP59BN3whVJ61XIkST).

 * Lecture slides will be made available online

 * [Exercises](https://microbiome.github.io/OMA/exercises.html)

 * [Resources](https://microbiome.github.io/OMA/resources.html)


## Installing and loading the required R packages {#packages}

This section shows how to install and load all required packages into
the R session, if needed. 


```r
# List of packages that we need from cran and bioc 
cran_pkg <- c("BiocManager", "bookdown", "dplyr", "ecodist", "ggplot2", 
              "gridExtra", "kableExtra",  "knitr", "scales", "vegan", "matrixStats")
bioc_pkg <- c("yulab.utils","ggtree","ANCOMBC", "ape", "DESeq2", "DirichletMultinomial", "mia", "miaViz")

# Get those packages that are already installed
cran_pkg_already_installed <- cran_pkg[ cran_pkg %in% installed.packages() ]
bioc_pkg_already_installed <- bioc_pkg[ bioc_pkg %in% installed.packages() ]

# Get those packages that need to be installed
cran_pkg_to_be_installed <- setdiff(cran_pkg, cran_pkg_already_installed)
bioc_pkg_to_be_installed <- setdiff(bioc_pkg, bioc_pkg_already_installed)

# Reorders bioc packages, so that mia and miaViz are first
bioc_pkg <- c(bioc_pkg[ bioc_pkg %in% c("mia", "miaViz") ], 
              bioc_pkg[ !bioc_pkg %in% c("mia", "miaViz") ] ) 

# Combine to one vector
packages <- c(bioc_pkg, cran_pkg)
packages_to_install <- c( bioc_pkg_to_be_installed, cran_pkg_to_be_installed )
```


```r
# If there are packages that need to be installed, install them 
if( length(packages_to_install) ) {
   BiocManager::install(packages_to_install)
}
```

Now all required packages are installed, so let's load them into the
session.  Some function names occur in multiple packages. That is why
we have prioritized the packages mia and miaViz on the list. Packages
that are loaded first have a higher priority.


```r
# Loading all packages into session. Returns true if package was successfully loaded.
loaded <- sapply(packages, require, character.only = TRUE)
as.data.frame(loaded)
```

```
##                      loaded
## mia                    TRUE
## miaViz                 TRUE
## yulab.utils            TRUE
## ggtree                 TRUE
## ANCOMBC                TRUE
## ape                    TRUE
## DESeq2                 TRUE
## DirichletMultinomial   TRUE
## BiocManager            TRUE
## bookdown               TRUE
## dplyr                  TRUE
## ecodist                TRUE
## ggplot2                TRUE
## gridExtra              TRUE
## kableExtra             TRUE
## knitr                  TRUE
## scales                 TRUE
## vegan                  TRUE
## matrixStats            TRUE
```

