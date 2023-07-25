# Microbiome Data {#containers}


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


## Data science framework

The building blocks of the framework are **data container**
(SummarizedExperiment and its derivatives), **packages** from various
developers using the TreeSE container, open **demonstration data
sets**, in a separate chapter \@ref(example-data), and **online
tutorials** including this online book as well as the various package
vignettes and other materials.


\includegraphics[width=18.67in]{general/figures/FigureOverviewV2_mod} 


## Data containers

`SummarizedExperiment` (`SE`) [@R_SummarizedExperiment] is a generic and highly optimized container for complex data
structures. It has become a common choice for analysing various types
of biomedical profiling data, such as RNAseq, ChIp-Seq, microarrays,
flow cytometry, proteomics, and single-cell
sequencing.

[`TreeSummarizedExperiment`] (`TreeSE`) [@R_TreeSummarizedExperiment] was developed as an extension to incorporate hierarchical
information (such as phylogenetic trees and sample hierarchies) and
reference sequences.

[`MultiAssayExperiment`] (`MAE`) [@Ramos2017] provides an organized way to bind several different data
structures together in a single object. For example, we can bind
microbiome data (in `TreeSE` format) with metabolomic profiling data
(in `SE`) format, with shared sample metadata. This is convenient and
robust for instance in subsetting and other data manipulation
tasks. Microbiome data can be part of multiomics experiments and
analysis strategies and we want to outline the understanding in which
we think the packages explained and used in this book relate to these
experiment layouts using the `TreeSummarizedExperiment` and classes
beyond.

This section provides an introductions to these data containers. In
microbiome data science, these containers link taxonomic abundance
tables with rich side information on the features and
samples. Taxonomic abundance data can be obtained by 16S rRNA amplicon
or metagenomic sequencing, phylogenetic microarrays, or by other
means. Many microbiome experiments include multiple versions and types
of data generated independently or derived from each other through
transformation or agglomeration. We start by providing recommendations
on how to represent different varieties of multi-table data within the
`TreeSummarizedExperiment` class.

The options and recommendations are summarized in Table \@ref(tab:options).


### Assay data {#assay-slot}

The original count-based taxonomic abundance tables may have different 
transformations, such as logarithmic, Centered Log-Ratio (CLR), or relative 
abundance. These are typically stored in _**assays**_.

Let us load example data and rename it as tse.


```r
library(mia)
data(hitchip1006, package="miaTime")
tse <- hitchip1006
```

The `assays` slot contains the experimental data as multiple count matrices. The result of `assays` is a list of matrices.


```r
assays(tse)
```

```
## List of length 1
## names(1): counts
```

Individual assays can be accessed via `assay`


```r
assay(tse, "counts")[1:5,1:7]
```

```
##                              Sample-1 Sample-2 Sample-3 Sample-4 Sample-5
## Actinomycetaceae                    0        0        0        0        0
## Aerococcus                          0        0        0        0        0
## Aeromonas                           0        0        0        0        0
## Akkermansia                        21       36      475       61       34
## Alcaligenes faecalis et rel.        1        1        1        2        1
##                              Sample-6 Sample-7
## Actinomycetaceae                    0        0
## Aerococcus                          0        0
## Aeromonas                           0        0
## Akkermansia                        14       27
## Alcaligenes faecalis et rel.        1        1
```

To illustrate the use of multiple assays, the relative abundance data can be 
calculated and stored along the original count data using `transformAssay`.


```r
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
assays(tse)
```

```
## List of length 2
## names(2): counts relabundance
```

Now there are two assays available in the `tse` object, `counts` and 
`relabundance`.


```r
assay(tse, "relabundance")[1:5,1:7]
```

```
##                               Sample-1  Sample-2  Sample-3  Sample-4  Sample-5
## Actinomycetaceae             0.0000000 0.000e+00 0.0000000 0.0000000 0.000e+00
## Aerococcus                   0.0000000 0.000e+00 0.0000000 0.0000000 0.000e+00
## Aeromonas                    0.0000000 0.000e+00 0.0000000 0.0000000 0.000e+00
## Akkermansia                  0.0027657 3.547e-03 0.0666106 0.0056195 2.833e-03
## Alcaligenes faecalis et rel. 0.0001317 9.854e-05 0.0001402 0.0001842 8.333e-05
##                               Sample-6  Sample-7
## Actinomycetaceae             0.0000000 0.0000000
## Aerococcus                   0.0000000 0.0000000
## Aeromonas                    0.0000000 0.0000000
## Akkermansia                  0.0017690 0.0045570
## Alcaligenes faecalis et rel. 0.0001264 0.0001688
```


Here the dimension of the count data remains unchanged in
transformation. This is in fact, a requirement for the assays.



### colData

`colData` contains data on the samples.


```r
colData(tse)
```

```
## DataFrame with 1151 rows and 10 columns
##                   age      sex nationality DNA_extraction_method  project
##             <integer> <factor>    <factor>              <factor> <factor>
## Sample-1           28   male            US                    NA        1
## Sample-2           24   female          US                    NA        1
## Sample-3           52   male            US                    NA        1
## Sample-4           22   female          US                    NA        1
## Sample-5           25   female          US                    NA        1
## ...               ...      ...         ...                   ...      ...
## Sample-1168        50   female Scandinavia                     r       40
## Sample-1169        31   female Scandinavia                     r       40
## Sample-1170        31   female Scandinavia                     r       40
## Sample-1171        52   male   Scandinavia                     r       40
## Sample-1172        52   male   Scandinavia                     r       40
##             diversity   bmi_group  subject      time      sample
##             <numeric>    <factor> <factor> <numeric> <character>
## Sample-1         5.76 severeobese        1         0    Sample-1
## Sample-2         6.06 obese              2         0    Sample-2
## Sample-3         5.50 lean               3         0    Sample-3
## Sample-4         5.87 underweight        4         0    Sample-4
## Sample-5         5.89 lean               5         0    Sample-5
## ...               ...         ...      ...       ...         ...
## Sample-1168      5.87 severeobese      244       8.1 Sample-1168
## Sample-1169      5.87 overweight       245       2.3 Sample-1169
## Sample-1170      5.92 overweight       245       8.2 Sample-1170
## Sample-1171      6.04 overweight       246       2.1 Sample-1171
## Sample-1172      5.74 overweight       246       7.9 Sample-1172
```



### rowData

`rowData` contains data on the features of the analyzed samples. Of particular
interest to the microbiome field, this is used to store taxonomic information.


```r
rowData(tse)
```

```
## DataFrame with 130 rows and 3 columns
##                                       Phylum          Family
##                                  <character>     <character>
## Actinomycetaceae              Actinobacteria  Actinobacteria
## Aerococcus                        Firmicutes         Bacilli
## Aeromonas                     Proteobacteria  Proteobacteria
## Akkermansia                  Verrucomicrobia Verrucomicrobia
## Alcaligenes faecalis et rel.  Proteobacteria  Proteobacteria
## ...                                      ...             ...
## Vibrio                        Proteobacteria  Proteobacteria
## Weissella et rel.                 Firmicutes         Bacilli
## Wissella et rel.                  Firmicutes         Bacilli
## Xanthomonadaceae              Proteobacteria  Proteobacteria
## Yersinia et rel.              Proteobacteria  Proteobacteria
##                                               Genus
##                                         <character>
## Actinomycetaceae                   Actinomycetaceae
## Aerococcus                               Aerococcus
## Aeromonas                                 Aeromonas
## Akkermansia                             Akkermansia
## Alcaligenes faecalis et rel. Alcaligenes faecalis..
## ...                                             ...
## Vibrio                                       Vibrio
## Weissella et rel.                 Weissella et rel.
## Wissella et rel.                   Wissella et rel.
## Xanthomonadaceae                   Xanthomonadaceae
## Yersinia et rel.                   Yersinia et rel.
```

### rowTree  

Phylogenetic trees also play an important role in the microbiome field. The 
`TreeSummarizedExperiment` class can keep track of features and node
relations via two functions, `rowTree` and `rowLinks`.

A tree can be accessed via `rowTree` as `phylo` object.       

```r
rowTree(tse)
```

```
## NULL
```

The links to the individual features are available through `rowLinks`.


```r
rowLinks(tse)
```

```
## NULL
```

Please note that there can be a 1:1 relationship between tree nodes and 
features, but this is not a must-have. This means there can be features, which
are not linked to nodes, and nodes, which are not linked to features. To change
the links in an existing object, the `changeTree` function is available.



### Alternative experiments {#alt-exp}

_**Alternative experiments**_ complement _assays_. They can contain
complementary data, which is no longer tied to the same dimensions as
the assay data. However, the number of samples (columns) must be the
same.

This can come into play, for instance, when one has taxonomic
abundance profiles quantified with different measurement technologies,
such as phylogenetic microarrays, amplicon sequencing, or metagenomic
sequencing. Another common use case is including abundance tables for
different taxonomic ranks. Such alternative experiments concerning the
same set of samples can be stored as

1. Separate _assays_ assuming that the taxonomic information can be mapped 
between features directly 1:1; or 
2. Data in the _altExp_ slot of the `TreeSummarizedExperiment`, if the feature 
dimensions differ. Each element of the _altExp_ slot is a `SummarizedExperiment`
or an object from a derived class with independent feature data.

The following shows how to store taxonomic abundance tables
agglomerated at different taxonomic levels. However, the data could as
well originate from entirely different measurement sources as long as
the samples match.

Let us first agglomerate the data to Phylum level. This yields a new
TreeSE data object.



```r
tse_phylum <- agglomerateByRank(tse, "Phylum", na.rm=TRUE)
# Both have the same number of columns (samples)
dim(tse)
```

```
## [1]  130 1151
```

```r
dim(tse_phylum)
```

```
## [1]    8 1151
```


Then we can add the new phylum-level data object as an alternative experiment in the original data.


```r
# Add the new data object to the original data object as an alternative experiment with the name "Phylum"
altExp(tse, "Phylum") <- tse_phylum

# Check the alternative experiment names available in the data
altExpNames(tse)
```

```
## [1] "Phylum"
```

We can now subset the data, for instance, and this acts on both altExp and assay data.


```r
tse[,1:10]
```

```
## class: TreeSummarizedExperiment 
## dim: 130 10 
## metadata(0):
## assays(2): counts relabundance
## rownames(130): Actinomycetaceae Aerococcus ... Xanthomonadaceae
##   Yersinia et rel.
## rowData names(3): Phylum Family Genus
## colnames(10): Sample-1 Sample-2 ... Sample-9 Sample-10
## colData names(10): age sex ... time sample
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(1): Phylum
## rowLinks: NULL
## rowTree: NULL
## colLinks: NULL
## colTree: NULL
```

```r
dim(altExp(tse[,1:10],"Phylum"))
```

```
## [1]  8 10
```

For more details on _altExp_, you can check the [introduction](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) to the `SingleCellExperiment` package [@R_SingleCellExperiment].



### MultiAssayExperiments {#mae}

_**Multiple experiments**_ relate to complementary measurement types,
such as transcriptomic or metabolomic profiling of the microbiome or
the host. Multiple experiments can be represented using the same
options as alternative experiments, or by using the
`MultiAssayExperiment` class [@Ramos2017]. Depending on how the 
datasets relate to each other the data can be stored as:

1. Separate _altExp_ if the samples can be matched directly 1:1; or
2. As `MultiAssayExperiment` objects, in which the connections between
samples are defined through a `sampleMap`. Each element on the
`experimentsList` of an `MultiAssayExperiment` is `matrix` or
`matrix`-like objects, including `SummarizedExperiment` objects, and 
the number of samples can differ between the elements.


<!--

```r
#TODO: Find the right dataset to explain a non 1:1 sample relationship
```
-->

For information have a look at the [intro vignette](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html) of the `MultiAssayExperiment` package.  

 
   Option   Rows (features)    Cols (samples)               Recommended  
---------   --------------    ---------------  ------------------------
   assays  	     match              match       Data transformations  
   altExp             free              match    Alternative experiments  
MultiAssay            free      free (mapping)    Multi-omic experiments    

Table: (\#tab:options) **Recommended options for storing multiple data tables in microbiome studies** The _assays_ are best suited for data transformations (one-to-one match between samples and columns across the assays). The _alternative experiments_ are particularly suitable for alternative versions of the data that are of same type but may have a different number of features (e.g. taxonomic groups); this is for instance the case with taxonomic abundance tables agglomerated at different levels (e.g. genus vs. phyla) or alternative profiling technologies (e.g. amplicon sequencing vs. shallow shotgun metagenomics). For alternative experiments one-to-one match between samples (cols) is libraryd but the alternative experiment tables can have different numbers of features (rows). Finally, elements of the _MultiAssayExperiment_ provide the most flexible way to incorporate multi-omic data tables with flexible numbers of samples and features. We recommend these conventions as the basis for methods development and application in microbiome studies.




## Demonstration data {#example-data}

Open demonstration data for testing and benchmarking purposes is
available from multiple locations. This chapter introduces some
options. The other chapters of this book provide ample examples about
the use of the data.




### Package data {#package-data}

The `mia` R package contains example datasets that are direct
conversions from the alternative `phyloseq` container to the
`TreeSummarizedExperiment` container.

List the [available
datasets](https://microbiome.github.io/mia/reference/index.html) in
the `mia` package:



```r
library(mia)
data(package="mia")
```

Load the `GlobalPatterns` data from the `mia` package:


```r
data("GlobalPatterns", package="mia")
GlobalPatterns
```

```
## class: TreeSummarizedExperiment 
## dim: 19216 26 
## metadata(0):
## assays(1): counts
## rownames(19216): 549322 522457 ... 200359 271582
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(26): CL3 CC1 ... Even2 Even3
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (19216 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```


#### HintikkaXOData {#hintikka-desc}

[HintikkaXOData](https://microbiome.github.io/microbiomeDataSets/reference/HintikkaXOData.html)
is derived from a study about the effects of fat diet and prebiotics on the
microbiome of rat models [@Hintikka2021]. It is available in the MAE data
container for R. The dataset is briefly presented in
[these slides](https://microbiome.github.io/outreach/hintikkaxo_presentation.html).


### ExperimentHub data

[ExperimentHub](https://bioconductor.org/packages/release/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html)
provides a variety of data resources, including the
[microbiomeDataSets](https://bioconductor.org/packages/release/data/experiment/html/microbiomeDataSets.html)
package [@Morgan2021; @microlahti2021].

A table of the available datasets is available through the
`availableDataSets` function.


```r
library(microbiomeDataSets)
availableDataSets()
```

```
##             Dataset
## 1  GrieneisenTSData
## 2    HintikkaXOData
## 3       LahtiMLData
## 4        LahtiMData
## 5       LahtiWAData
## 6      OKeefeDSData
## 7 SilvermanAGutData
## 8        SongQAData
## 9   SprockettTHData
```

All data are downloaded from ExperimentHub and cached for local
re-use. Check the [man pages of each
function](https://microbiome.github.io/microbiomeDataSets/reference/index.html)
for a detailed documentation of the data contents and references. Let
us retrieve a *[MultiAssayExperiment](https://bioconductor.org/packages/3.17/MultiAssayExperiment)* dataset:


```r
# mae <- HintikkaXOData()
# Since HintikkaXOData is now added to mia, we can load it directly from there
# We suggest to check other datasets from microbiomeDataSets
data(HintikkaXOData, package = "mia")
mae <- HintikkaXOData
```

Data is available in *[SummarizedExperiment](https://bioconductor.org/packages/3.17/SummarizedExperiment)*, `r
Biocpkg("TreeSummarizedExperiment")` and `r
Biocpkg("MultiAssayExperiment")` data containers; see the separate
page on [alternative
containers](https://microbiome.github.io/OMA/multitable.html) for more
details.


### Curated metagenomic data

[curatedMetagenomicData](https://bioconductor.org/packages/release/data/experiment/html/curatedMetagenomicData.html)
is a large collection of curated human microbiome datasets, provided as
`(Tree)SummarizedExperiment` objects [@Pasolli2017]. The resource
provides curated human microbiome data including gene families, marker
abundance, marker presence, pathway abundance, pathway coverage, and
relative abundance for samples from different body sites. See the
package homepage for more details on data availability and access.

As one example, let us retrieve the Vatanen (2016) [@Vatanen2016] data
set. This is a larger collection with a bit longer download time.


```r
library(curatedMetagenomicData)
tse <- curatedMetagenomicData("Vatanen*", dryrun = FALSE, counts = TRUE)
```



### Other data sources

The current collections provide access to vast microbiome data
resources. The output has to be converted into TreeSE/MAE separately.

- [MGnifyR](https://github.com/beadyallen/MGnifyR) provides access to [EBI/MGnify](https://www.ebi.ac.uk/metagenomics/) 
- [qiitr](https://github.com/cran/qiitr) provides access to [QIITA](https://qiita.com/about) 


## Loading experimental microbiome data

### 16S workflow

Result of amplicon sequencing is a large number of files that include all the sequences
that were read from samples. Those sequences need to be matched with taxa. Additionally,
we need to know how many times each taxa were found from each sample. 

There are several algorithms to do that, and DADA2 is one of the most common. 
You can find DADA2 pipeline tutorial, for example, 
[here](https://benjjneb.github.io/dada2/tutorial.html).
After the DADA2 portion of the tutorial is completed, the data is stored into _phyloseq_ object 
(Bonus: Handoff to phyloseq). To store the data to _TreeSummarizedExperiment_,
follow the example below. 

You can find full workflow script without further explanations and comments from 
[here](https://github.com/microbiome/OMA/blob/master/dada2_workflow.Rmd)



Load required packages.


```r
library(mia)
library(ggplot2)
library(BiocManager)
library(Biostrings)
```

Create arbitrary example sample metadata like it was done in the tutorial. Usually, 
sample metadata is imported as a file.


```r
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
```

Convert data into right format and create a _TreeSE_ object.


```r
# Create a list that contains assays
counts <- t(seqtab.nochim)
counts <- as.matrix(counts)
assays <- SimpleList(counts = counts)

# Convert colData and rowData into DataFrame
samdf <- DataFrame(samdf)
taxa <- DataFrame(taxa)

# Create TreeSE
tse <- TreeSummarizedExperiment(assays = assays,
                                colData = samdf,
                                rowData = taxa
                                )

# Remove mock sample like it is also done in DADA2 pipeline tutorial
tse <- tse[ , colnames(tse) != "mock"]
```

Add sequences into _referenceSeq_ slot and convert rownames into simpler format.


```r
# Convert sequences into right format
dna <- Biostrings::DNAStringSet( rownames(tse) )
# Add sequences into referenceSeq slot
referenceSeq(tse) <- dna
# Convert rownames into ASV_number format
rownames(tse) <- paste0("ASV", seq( nrow(tse) ))
tse
```

```
## class: TreeSummarizedExperiment 
## dim: 232 20 
## metadata(0):
## assays(1): counts
## rownames(232): ASV1 ASV2 ... ASV231 ASV232
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(20): F3D0 F3D1 ... F3D9 Mock
## colData names(4): Subject Gender Day When
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: NULL
## rowTree: NULL
## colLinks: NULL
## colTree: NULL
## referenceSeq: a DNAStringSet (232 sequences)
```



### Import from external files {#import-from-file}

Microbiome (taxonomic) profiling data is commonly distributed in
various file formats. You can import such external data files as a
(Tree)SummarizedExperiment object, but the details depend on the file
format. Here, we provide examples for common formats. Some datasets and raw
files to learn how to import raw data and construct TreeSE/MAE containers are
available in the
[microbiome data repository](https://github.com/microbiome/data).


#### CSV import

**CSV data tables** can be imported with the standard R functions,
  then converted to the desired format. For detailed examples, you can
  check the [Bioconductor course
  material](https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html)
  by Martin Morgan. You can also check the [example
  files](https://github.com/microbiome/OMA/tree/master/data) and
  construct your own CSV files accordingly.

Recommendations for the CSV files are the following. File names are
arbitrary; we refer here to the same names as in the examples:

- Abundance table (`assay_taxa.csv`): data matrix (features x
  samples); first column provides feature IDs, the first row provides
  sample IDs; other values should be numeric (abundances).

- Row data (`rowdata_taxa.csv`): data table (features x info); first
  column provides feature IDs, the first row provides column headers;
  this file usually contains the taxonomic mapping between different
  taxonomic levels. Ideally, the feature IDs (row names) match one-to-one with
  the abundance table row names. 

- Column data (`coldata.csv`): data table (samples x info); first
  column provides sample IDs, the first row provides column headers;
  this file usually contains the sample metadata/phenodata (such as
  subject age, health etc). Ideally, the sample IDs match one-to-one with
  the abundance table column names. 

After you have set up the CSV files, you can read them in R:


```r
count_file  <- "data/assay_taxa.csv"
tax_file    <- "data/rowdata_taxa.csv"
sample_file <- "data/coldata.csv"

# Load files
counts  <- read.csv(count_file, row.names=1)   # Abundance table (e.g. ASV data; to assay data)
tax     <- read.csv(tax_file, row.names=1)     # Taxonomy table (to rowData)
samples <- read.csv(sample_file, row.names=1)  # Sample data (to colData)
```

After reading the data in R, ensure the following:

- abundance table (`counts`): numeric `matrix`, with feature IDs as
  rownames and sample IDs as column names

- rowdata (`tax`): `DataFrame`, with feature IDs as rownames. If this
  is a `data.frame` you can use the function `DataFrame()` to change
  the format. Column names are free but in microbiome analysis they
  usually they refer to taxonomic ranks. The rownames in rowdata
  should match with rownames in abundance table.

- coldata (`samples`): `DataFrame`, with sample IDs as rownames. If
  this is a `data.frame` you can use the function `DataFrame()` to
  change the format.  Column names are free. The rownames in coldata
  should match with colnames in abundance table.

**Always ensure that the tables have rownames!** The _TreeSE_ constructor compares 
rownames and ensures that, for example, right samples are linked with right patient.

Also ensure that the row and column names match one-to-one between
abundance table, rowdata, and coldata:


```r
# Match rows and columns
counts <- counts[rownames(tax), rownames(samples)]

# Let us ensure that the data is in correct (numeric matrix) format:
counts <- as.matrix(counts)
```

If you hesitate about the format of the data, you can compare to one
of the available demonstration datasets, and make sure that your data
components have the same format.

There are many different source files and many different ways to read
data in R. One can do data manipulation in R as well. Investigate the
entries as follows.



```r
# coldata rownames match assay colnames
all(rownames(samples) == colnames(counts)) # our dataset
```

```
## [1] TRUE
```

```r
class(samples) # should be data.frame or DataFrame
```

```
## [1] "data.frame"
```

```r
# rowdata rownames match assay rownames
all(rownames(tax) == rownames(counts)) # our dataset
```

```
## [1] TRUE
```

```r
class(tax) # should be data.frame or DataFrame
```

```
## [1] "data.frame"
```

```r
# Counts 
class(counts) # should be a numeric matrix
```

```
## [1] "matrix" "array"
```


### Constructing TreeSummarizedExperiment

Now let us create the TreeSE object from the input data tables. Here
we also convert the data objects in their preferred formats:

   - counts --> numeric matrix
   - rowData --> DataFrame
   - colData --> DataFrame

The `SimpleList` could be used to include multiple alternative assays, if
necessary.



```r
# Create a TreeSE
tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = counts),
                                     colData = DataFrame(samples),
                                     rowData = DataFrame(tax))

tse_taxa
```

```
## class: TreeSummarizedExperiment 
## dim: 12706 40 
## metadata(0):
## assays(1): counts
## rownames(12706): GAYR01026362.62.2014 CVJT01000011.50.2173 ...
##   JRJTB:03787:02429 JRJTB:03787:02478
## rowData names(7): Phylum Class ... Species OTU
## colnames(40): C1 C2 ... C39 C40
## colData names(6): Sample Rat ... Fat XOS
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: NULL
## rowTree: NULL
## colLinks: NULL
## colTree: NULL
```

Now you should have a ready-made TreeSE data object that can be used in downstream analyses.


### Constructing MultiAssayExperiment

To construct a _MultiAssayExperiment_ object, just combine multiple _TreeSE_ data containers. 
Here we import metabolite data from the same study.


```r
count_file <- "data/assay_metabolites.csv"
sample_file <- "data/coldata.csv"

# Load files
counts  <- read.csv(count_file, row.names=1)  
samples <- read.csv(sample_file, row.names=1)

# Create a TreeSE for the metabolite data
tse_metabolite <- TreeSummarizedExperiment(assays = SimpleList(concs = as.matrix(counts)),
                                           colData = DataFrame(samples))

tse_metabolite
```

```
## class: TreeSummarizedExperiment 
## dim: 38 40 
## metadata(0):
## assays(1): concs
## rownames(38): Butyrate Acetate ... Malonate 1,3-dihydroxyacetone
## rowData names(0):
## colnames(40): C1 C2 ... C39 C40
## colData names(6): Sample Rat ... Fat XOS
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: NULL
## rowTree: NULL
## colLinks: NULL
## colTree: NULL
```

Now we can combine these two experiments into _MAE_.


```r
# Create an ExperimentList that includes experiments
experiments <- ExperimentList(microbiome = tse_taxa, 
                              metabolite = tse_metabolite)

# Create a MAE
mae <- MultiAssayExperiment(experiments = experiments)

mae
```

```
## A MultiAssayExperiment object of 2 listed
##  experiments with user-defined names and respective classes.
##  Containing an ExperimentList class object of length 2:
##  [1] microbiome: TreeSummarizedExperiment with 12706 rows and 40 columns
##  [2] metabolite: TreeSummarizedExperiment with 38 rows and 40 columns
## Functionality:
##  experiments() - obtain the ExperimentList instance
##  colData() - the primary/phenotype DataFrame
##  sampleMap() - the sample coordination DataFrame
##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
##  *Format() - convert into a long or wide DataFrame
##  assays() - convert ExperimentList to a SimpleList of matrices
##  exportClass() - save data to flat files
```



### Import functions for standard formats

Specific import functions are provided for:

-   Biom files (see `help(mia::loadFromBiom)`)
-   QIIME2 files (see `help(mia::loadFromQIIME2)`)
-   Mothur files (see `help(mia::loadFromMothur)`)


#### Biom import

This example shows how [Biom files](https://biom-format.org/) are imported into
a `TreeSummarizedExperiment` object.

The data is from following publication: 
Tengeler AC _et al._ (2020) [**Gut microbiota from persons with
attention-deficit/hyperactivity disorder affects the brain in
mice**](https://doi.org/10.1186/s40168-020-00816-x). 

The dataset consists of 3 files:

-   biom file: abundance table and taxonomy information
-   csv file: sample metadata
-   tree file: phylogenetic tree


Store the data in your desired local directory (for instance, _data/_ under the
working directory), and define source file paths


```r
biom_file_path <- "data/Aggregated_humanization2.biom"
sample_meta_file_path <- "data/Mapping_file_ADHD_aggregated.csv"
tree_file_path <- "data/Data_humanization_phylo_aggregation.tre"
```

Now we can load the biom data into a SummarizedExperiment (SE) object.


```r
library(mia)

# Imports the data
se <- loadFromBiom(biom_file_path)

# Check
se
```

```
## class: TreeSummarizedExperiment 
## dim: 151 27 
## metadata(0):
## assays(1): counts
## rownames(151): 1726470 1726471 ... 17264756 17264757
## rowData names(6): taxonomy1 taxonomy2 ... taxonomy5 taxonomy6
## colnames(27): A110 A111 ... A38 A39
## colData names(0):
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: NULL
## rowTree: NULL
## colLinks: NULL
## colTree: NULL
```

The `assays` slot includes a list of abundance tables. The imported
abundance table is named as "counts".  Let us inspect only the first
cols and rows.


```r
assays(se)$counts[1:3, 1:3]
```

```
##           A110  A111  A12
## 1726470  17722 11630    0
## 1726471  12052     0 2679
## 17264731     0   970    0
```

The `rowdata` includes taxonomic information from the biom file. The `head()` command
shows just the beginning of the data table for an overview.

`knitr::kable()` is for printing the information more nicely.


```r
head(rowData(se))
```

```
## DataFrame with 6 rows and 6 columns
##             taxonomy1          taxonomy2           taxonomy3
##           <character>        <character>         <character>
## 1726470  "k__Bacteria   p__Bacteroidetes      c__Bacteroidia
## 1726471  "k__Bacteria   p__Bacteroidetes      c__Bacteroidia
## 17264731 "k__Bacteria   p__Bacteroidetes      c__Bacteroidia
## 17264726 "k__Bacteria   p__Bacteroidetes      c__Bacteroidia
## 1726472  "k__Bacteria p__Verrucomicrobia c__Verrucomicrobiae
## 17264724 "k__Bacteria   p__Bacteroidetes      c__Bacteroidia
##                      taxonomy4              taxonomy5           taxonomy6
##                    <character>            <character>         <character>
## 1726470       o__Bacteroidales      f__Bacteroidaceae     g__Bacteroides"
## 1726471       o__Bacteroidales      f__Bacteroidaceae     g__Bacteroides"
## 17264731      o__Bacteroidales  f__Porphyromonadaceae g__Parabacteroides"
## 17264726      o__Bacteroidales      f__Bacteroidaceae     g__Bacteroides"
## 1726472  o__Verrucomicrobiales f__Verrucomicrobiaceae     g__Akkermansia"
## 17264724      o__Bacteroidales      f__Bacteroidaceae     g__Bacteroides"
```

These taxonomic rank names (column names) are not real rank
names. Let’s replace them with real rank names.

In addition to that, the taxa names include, e.g., '"k__' before the name, so let's
make them cleaner by removing them. 


```r
names(rowData(se)) <- c("Kingdom", "Phylum", "Class", "Order", 
                        "Family", "Genus")

# Goes through the whole DataFrame. Removes '.*[kpcofg]__' from strings, where [kpcofg] 
# is any character from listed ones, and .* any character.
rowdata_modified <- BiocParallel::bplapply(rowData(se), 
                                           FUN = stringr::str_remove, 
                                           pattern = '.*[kpcofg]__')

# Genus level has additional '\"', so let's delete that also
rowdata_modified <- BiocParallel::bplapply(rowdata_modified, 
                                           FUN = stringr::str_remove, 
                                           pattern = '\"')

# rowdata_modified is a list, so it is converted back to DataFrame format. 
rowdata_modified <- DataFrame(rowdata_modified)

# And then assigned back to the SE object
rowData(se) <- rowdata_modified

# Now we have a nicer table
head(rowData(se))
```

```
## DataFrame with 6 rows and 6 columns
##              Kingdom          Phylum            Class              Order
##          <character>     <character>      <character>        <character>
## 1726470     Bacteria   Bacteroidetes      Bacteroidia      Bacteroidales
## 1726471     Bacteria   Bacteroidetes      Bacteroidia      Bacteroidales
## 17264731    Bacteria   Bacteroidetes      Bacteroidia      Bacteroidales
## 17264726    Bacteria   Bacteroidetes      Bacteroidia      Bacteroidales
## 1726472     Bacteria Verrucomicrobia Verrucomicrobiae Verrucomicrobiales
## 17264724    Bacteria   Bacteroidetes      Bacteroidia      Bacteroidales
##                       Family           Genus
##                  <character>     <character>
## 1726470       Bacteroidaceae     Bacteroides
## 1726471       Bacteroidaceae     Bacteroides
## 17264731  Porphyromonadaceae Parabacteroides
## 17264726      Bacteroidaceae     Bacteroides
## 1726472  Verrucomicrobiaceae     Akkermansia
## 17264724      Bacteroidaceae     Bacteroides
```

We notice that the imported biom file did not contain the sample meta data
yet, so it includes an empty data frame.


```r
head(colData(se))
```

```
## DataFrame with 6 rows and 0 columns
```

Let us add a sample metadata file.


```r
# We use this to check what type of data it is
# read.table(sample_meta_file_path)

# It seems like a comma separated file and it does not include headers
# Let us read it and then convert from data.frame to DataFrame
# (required for our purposes)
sample_meta <- DataFrame(read.table(sample_meta_file_path, sep = ",", header = FALSE))

# Add sample names to rownames
rownames(sample_meta) <- sample_meta[,1]

# Delete column that included sample names
sample_meta[,1] <- NULL

# We can add headers
colnames(sample_meta) <- c("patient_status", "cohort", "patient_status_vs_cohort", "sample_name")

# Then it can be added to colData
colData(se) <- sample_meta
```

Now `colData` includes the sample metadata.


```r
head(colData(se))
```

```
## DataFrame with 6 rows and 4 columns
##      patient_status      cohort patient_status_vs_cohort sample_name
##         <character> <character>              <character> <character>
## A110           ADHD    Cohort_1            ADHD_Cohort_1        A110
## A12            ADHD    Cohort_1            ADHD_Cohort_1         A12
## A15            ADHD    Cohort_1            ADHD_Cohort_1         A15
## A19            ADHD    Cohort_1            ADHD_Cohort_1         A19
## A21            ADHD    Cohort_2            ADHD_Cohort_2         A21
## A23            ADHD    Cohort_2            ADHD_Cohort_2         A23
```

Now, let's add a phylogenetic tree.

The current data object, se, is a SummarizedExperiment object. This
does not include a slot for adding a phylogenetic tree. In order to do
this, we can convert the SE object to an extended TreeSummarizedExperiment
object which includes also a `rowTree` slot.

TreeSummarizedExperiment contains also other additional slots and features which
is why we recommend to use `TreeSE`.


```r
tse <- as(se, "TreeSummarizedExperiment")

# tse includes same data as se
tse
```

```
## class: TreeSummarizedExperiment 
## dim: 151 27 
## metadata(0):
## assays(1): counts
## rownames(151): 1726470 1726471 ... 17264756 17264757
## rowData names(6): Kingdom Phylum ... Family Genus
## colnames(27): A110 A12 ... A35 A38
## colData names(4): patient_status cohort patient_status_vs_cohort
##   sample_name
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: NULL
## rowTree: NULL
## colLinks: NULL
## colTree: NULL
```

Next, let us read the tree data file and add it to the R data object (tse).


```r
# Reads the tree file
tree <- ape::read.tree(tree_file_path)

# Add tree to rowTree
rowTree(tse) <- tree

# Check
tse
```

```
## class: TreeSummarizedExperiment 
## dim: 151 27 
## metadata(0):
## assays(1): counts
## rownames(151): 1726470 1726471 ... 17264756 17264757
## rowData names(6): Kingdom Phylum ... Family Genus
## colnames(27): A110 A12 ... A35 A38
## colData names(4): patient_status cohort patient_status_vs_cohort
##   sample_name
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (151 rows)
## rowTree: 1 phylo tree(s) (151 leaves)
## colLinks: NULL
## colTree: NULL
```

Now `rowTree` includes a phylogenetic tree:


```r
head(rowTree(tse))
```


### Conversions between data formats in R

If the data has already been imported in R in another format, it
can be readily converted into `TreeSummarizedExperiment`, as shown in our next
example. Note that similar conversion functions to
`TreeSummarizedExperiment` are available for multiple data formats via
the `mia` package (see makeTreeSummarizedExperimentFrom* for phyloseq,
Biom, and DADA2).


```r
library(mia)

# phyloseq example data
data(GlobalPatterns, package="phyloseq") 
GlobalPatterns_phyloseq <- GlobalPatterns
GlobalPatterns_phyloseq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
## sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
## tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
```


```r
# convert phyloseq to TSE
GlobalPatterns_TSE <- makeTreeSummarizedExperimentFromPhyloseq(GlobalPatterns_phyloseq) 
GlobalPatterns_TSE
```

```
## class: TreeSummarizedExperiment 
## dim: 19216 26 
## metadata(0):
## assays(1): counts
## rownames(19216): 549322 522457 ... 200359 271582
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(26): CL3 CC1 ... Even2 Even3
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (19216 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```

We can also convert `TreeSummarizedExperiment` objects into `phyloseq`
with respect to the shared components that are supported by both
formats (i.e. taxonomic abundance table, sample metadata, taxonomic
table, phylogenetic tree, sequence information). This is useful for
instance when additional methods are available for `phyloseq`.


```r
# convert TSE to phyloseq
GlobalPatterns_phyloseq2 <- makePhyloseqFromTreeSummarizedExperiment(GlobalPatterns_TSE) 
GlobalPatterns_phyloseq2
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
## sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
## tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
```


Conversion is possible between other data formats. Interested readers can refer to the following functions:
* [makeTreeSummarizedExperimentFromDADA2](https://microbiome.github.io/mia/reference/makeTreeSummarizedExperimentFromDADA2.html)
* [makeSummarizedExperimentFromBiom](https://microbiome.github.io/mia/reference/makeSummarizedExperimentFromBiom.html)
* [loadFromMetaphlan](https://microbiome.github.io/mia/reference/loadFromMetaphlan.html)
* [readQZA](https://microbiome.github.io/mia/reference/loadFromQIIME2.html)



