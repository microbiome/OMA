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


The Bioconductor microbiome data science framework consists of:

- **data containers**, designed to organize multi-assay microbiome data
- **R/Bioconductor packages** that provide dedicated methods 
- **community** of users and developers 

<img src="general/figures/ecosystem.png" width="100" alt="mia logo" align="right" style="margin: 0 1em 0 1em" />

This section provides an overview of the package ecosystem. Section
\@ref(example-data) links to various open microbiome data resources
that support this framework.


## Package installation

You can install all packages that are required to run every example in this book via the following command:


```r
source("https://raw.githubusercontent.com/microbiome/OMA/master/install_packages.R")
```

### Installing specific packages {#packages_specific}

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


## Package ecosystem {#ecosystem}

Methods for `(Tree)SummarizedExperiment` and `MultiAssayExperiment`
data containers are provided by multiple independent developers
through R/Bioconductor packages. Some of these are listed below (tips
on new packages are [welcome](https://microbiome.github.io)).


### mia package family

The mia package family provides general methods for microbiome data wrangling, analysis and visualization. 

- [mia](https://microbiome.github.io/mia/): Microbiome analysis tools [@R_mia]
- [miaViz](https://microbiome.github.io/miaViz/): Microbiome analysis specific visualization [@Ernst2022]
- [miaSim](https://microbiome.github.io/miaSim/): Microbiome data simulations [@Simsek2021]
- [miaTime](https://microbiome.github.io/miaTime/): Microbiome time series analysis [@Lahti2021]


### Differential abundance {#sub-diff-abund}

The following DA methods support `(Tree)SummarizedExperiment`.

- [ANCOMBC](https://bioconductor.org/packages/devel/bioc/html/ANCOMBC.html) for differential abundance analysis
- [benchdamic](https://bioconductor.org/packages/release/bioc/vignettes/benchdamic/inst/doc/intro.html) for benchmarking differential abundance methods
- [ALDEx2](https://www.bioconductor.org/packages/release/bioc/html/ALDEx2.html) for differential abundance analysis


### Other packages

- [philr](http://bioconductor.org/packages/devel/bioc/html/philr.html) (@Silverman2017) phylogeny-aware phILR transformation
- [MicrobiotaProcess](https://bioconductor.org/packages/release/bioc/html/MicrobiotaProcess.html) for "tidy" analysis of microbiome and other ecological data
- [Tools for Microbiome
  Analysis](https://microsud.github.io/Tools-Microbiome-Analysis/)
  site listed over 130 R packages for microbiome data science in
  2023. Many of these are not in Bioconductor, or do not directly
  support the data containers used in this book but can be often used
  with minor modifications. The phyloseq-based tools can be used by
  converting the TreeSE data into phyloseq with
  `makePhyloseqFromTreeSummarizedExperiment`.


### Open microbiome data 

Hundreds of published microbiome data sets are readily available in
these data containers (see \@ref(example-data)).

