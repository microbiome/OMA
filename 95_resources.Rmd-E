
# Resources {#resources}

## Data containers

### Resources for TreeSummarizedExperiment

 * SingleCellExperiment [@R-SingleCellExperiment]
   + [Online tutorial](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)
   + [Project page](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
 * SummarizedExperiment [@R-SummarizedExperiment]
   + [Online tutorial](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html)
   + [Project page](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
 * TreeSummarizedExperiment [@R-TreeSummarizedExperiment]
   + [Online tutorial](https://bioconductor.org/packages/release/bioc/vignettes/TreeSummarizedExperiment/inst/doc/Introduction_to_treeSummarizedExperiment.html)
   + [Project page](https://www.bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html)
   + Publication: [@Huang2021]
   

### Other relevant containers

* [DataFrame](https://rdrr.io/bioc/S4Vectors/man/DataFrame-class.html) which behaves similarly to `data.frame`, yet efficient and fast when used with large datasets.
* [DNAString](https://rdrr.io/bioc/Biostrings/man/DNAString-class.html) along with `DNAStringSet`,`RNAString` and `RNAStringSet`  efficient storage and handling of long biological sequences are offered within the Biostrings package [@R-Biostrings].
* GenomicRanges ([@GenomicRanges2013]) offers an efficient representation and manipulation of genomic annotations and alignments, see e.g. `GRanges` and `GRangesList` at [An Introduction to the GenomicRangesPackage](https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html).

[NGS Analysis Basics](http://girke.bioinformatics.ucr.edu/GEN242/tutorials/rsequences/rsequences/) provides a walk-through of the above-mentioned features with detailed examples.


### Alternative containers for microbiome data

The `phyloseq` package and class became the first widely used data
container for microbiome data science in R. Many methods for taxonomic
profiling data are readily available for this class. We provide here a
short description how `phyloseq` and `*Experiment` classes relate to
each other.

`assays` : This slot is similar to `otu_table` in `phyloseq`. In a
               `SummarizedExperiment` object multiple assays, raw
               counts, transformed counts can be stored. See also
               [@Ramos2017]
               for storing data from multiple `experiments` such as
               RNASeq, Proteomics, etc.  `rowData` : This slot is
               similar to `tax_table` in `phyloseq` to store taxonomic
               information.  `colData` : This slot is similar to
               `sample_data` in `phyloseq` to store information
               related to samples.  `rowTree` : This slot is similar
               to `phy_tree` in `phyloseq` to store phylogenetic tree.

In this book, you will come across terms like `FeatureIDs` and
`SampleIDs`.  `FeatureIDs` : These are basically OTU/ASV ids which are
row names in `assays` and `rowData`.  `SampleIDs` : As the name
suggests, these are sample ids which are column names in `assays` and
row names in `colData`.

`FeatureIDs` and `SampleIDs` are used but the technical terms
`rownames` and `colnames` are encouraged to be used, since they relate
to actual objects we work with.

<img
src="https://raw.githubusercontent.com/FelixErnst/TreeSummarizedExperiment
/2293440c6e70ae4d6e978b6fdf2c42fdea7fb36a/vignettes/tse2.png"
width="100%"/>


### Resources for phyloseq

The (Tree)SummarizedExperiment objects can be converted into the alternative phyloseq format, for which further methods are available.

 * [List of R tools for microbiome analysis](https://microsud.github.io/Tools-Microbiome-Analysis/)
 * phyloseq [@McMurdie2013]
 * [microbiome tutorial](http://microbiome.github.io/tutorials/)
 * [microbiomeutilities](https://microsud.github.io/microbiomeutilities/)
 * Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses [@Callahan2016].





## R programming resources

If you are new to R, you could try [swirl](https://swirlstats.com/)
for a kickstart to R programming. Further support resources are
available through the Bioconductor
project [@Huber2015].

 * R programming basics: [Base R](https://www.rstudio.com/wp-content/uploads/2016/10/r-cheat-sheet-3.pdf)
 * Basics of R programming: [Base R](https://raw.githubusercontent.com/rstudio/cheatsheets/master/base-r.pdf)
 * [R cheat sheets](https://www.rstudio.com/resources/cheatsheets/)
 * R visualization with [ggplot2](https://www.rstudio.com/wp-content/uploads/2016/11/ggplot2-cheatsheet-2.1.pdf) 
 * [R graphics cookbook](http://www.cookbook-r.com/Graphs/)

Rmarkdown

* Rmarkdown tips [@Xie2020]


RStudio

* [RStudio cheat sheet](https://www.rstudio.com/wp-content/uploads/2016/01/rstudio-IDE-cheatsheet.pdf) 



### Bioconductor Classes {#bioc_intro}

**S4 system**

S4 class system has brought several useful features to the
object-oriented programming paradigm within R, and it is constantly
deployed in R/Bioconductor packages [@Huber2015].
    
|   Online Document:

* Hervé Pagès, [A quick overview of the S4 class system](https://bioconductor.org/packages/release/bioc/vignettes/S4Vectors/inst/doc/S4QuickOverview.pdf).
* Laurent Gatto, [A practical tutorial on S4 programming](https://bioconductor.org/help/course-materials/2013/CSAMA2013/friday/afternoon/S4-tutorial.pdf)
* How S4 Methods Work [@Chambers2006]

|   Books:

* John M. Chambers. Software for Data Analysis: Programming with R. Springer, New York, 2008. ISBN-13 978-0387759357 [@Chambers2008]
* I Robert Gentleman. R Programming for Bioinformatics. Chapman & Hall/CRC, New York, 2008. ISBN-13 978-1420063677 [@gentleman2008r]



## Reproducible reporting with Quarto {#quarto}

Reproducible reporting is the starting point for robust interactive
data science. Perform the following tasks:

 * If you are entirely new to literate programming and reproducible
   reporting, take [this](https://www.markdowntutorial.com/) 10 minute
   tutorial to get introduced to the most important functions within
   Markdown. 

 * Create a Quarto template in RStudio, and render it into a
   document (markdown, PDF, docx or other format). In case you are new
   to Quarto, see [this page](https://quarto.org/).

 * Examples are tips for closely related Rmarkdown are available in
   the online tutorial to reproducible reporting by [Dr. C Titus
   Brown](https://rpubs.com/marschmi/RMarkdown); see also [Rmarkdown
   cheatsheet](https://www.rstudio.com/wp-content/uploads/2015/02/rmarkdown-cheatsheet.pdf)


**Figure sources:** 

**Original article**
-   Huang R _et al_. (2021) TreeSummarizedExperiment: a S4 class 
for data with hierarchical structure. F1000Research 9:1246. [@Huang2021]

**Reference Sequence slot extension**
- Lahti L _et al_. (2020) [Upgrading the R/Bioconductor ecosystem for microbiome 
research](https://doi.org/10.7490/
f1000research.1118447.1) F1000Research 9:1464 (slides).
