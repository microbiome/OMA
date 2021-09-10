# Introduction to Bioconductor Classes {-}

**S4 system**

S4 class system has brought several useful features to the
object-oriented programming paradigm within R.  It has been constantly
deployed in R packages, including the
[Bioconductor](https://bioconductor.org/) package.


|   Books:

* John M. Chambers. Software for Data Analysis: Programming with R. Springer, New York, 2008. ISBN-13 978-0387759357.
* I Robert Gentleman. R Programming for Bioinformatics. Chapman & Hall/CRC, New York, 2008. ISBN-13 978-1420063677
    
|   Online Document:

* Hervé Pagès, [A quick overview of the S4 class system](https://bioconductor.org/packages/release/bioc/vignettes/S4Vectors/inst/doc/S4QuickOverview.pdf).
* Laurent Gatto, [A practical tutorial on S4 programming](https://bioconductor.org/help/course-materials/2013/CSAMA2013/friday/afternoon/S4-tutorial.pdf)
* John M. Chambers. [How S4 Methods Work](http://developer.r-project.org/howMethodsWork.pdf)

**Additional useful classes include:**

* [DataFrame](https://rdrr.io/bioc/S4Vectors/man/DataFrame-class.html) which behaves similarly to `data.frame`, yet efficient and fast when used with large datasets.
* [DNAString](https://rdrr.io/bioc/Biostrings/man/DNAString-class.html) along with `DNAStringSet`,`RNAString` and `RNAStringSet`  efficient storage and handling of long biological sequences are offered within the [Biostrings](https://rdrr.io/bioc/Biostrings/) package.
* [GenomicRanges](https://bioconductor.org/packages/3.13/bioc/html/GenomicRanges.html) offers an efficient representation and manipulation of genomic annotations and alignments, see e.g. `GRanges` and `GRangesList` at [An Introduction to the GenomicRangesPackage](https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html).

[NGS Analysis Basics](http://girke.bioinformatics.ucr.edu/GEN242/tutorials/rsequences/rsequences/) is an example walk-through of the above mentioned features, with detailed examples provided.
