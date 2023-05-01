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
vignettes and other material.


\includegraphics[width=18.67in]{general/figures/FigureOverviewV2_mod} 


## Data containers

`SummarizedExperiment` (`SE`) [@R-SummarizedExperiment] is a generic and highly optimized container for complex data
structures. It has become a common choice for analysing various types
of biomedical profiling data, such as RNAseq, ChIp-Seq, microarrays,
flow cytometry, proteomics, and single-cell
sequencing.

[`TreeSummarizedExperiment`] (`TreeSE`) [@R-TreeSummarizedExperiment] was developed as an extension to incorporate hierarchical
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


### Assay data 

The original count-based taxonomic abundance tables may have different 
transformations, such as logarithmic, Centered Log-Ratio (CLR), or relative 
abundance. These are typically stored in _**assays**_.


```r
library(mia)
data(GlobalPatterns, package="mia")
tse <- GlobalPatterns
assays(tse)
```

```
## List of length 1
## names(1): counts
```

The `assays` slot contains the experimental data as count matrices. Multiple 
matrices can be stored the result of `assays` is actually a list of matrices.


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
##        CL3 CC1 SV1 M31Fcsw M11Fcsw M31Plmr M11Plmr
## 549322   0   0   0       0       0       0       0
## 522457   0   0   0       0       0       0       0
## 951      0   0   0       0       0       0       1
## 244423   0   0   0       0       0       0       0
## 586076   0   0   0       0       0       0       0
```

To illustrate the use of multiple assays, the relative abundance data can be 
calculated and stored along the original count data using `relAbundanceCounts`.


```r
tse <- relAbundanceCounts(tse)
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
##        CL3 CC1 SV1 M31Fcsw M11Fcsw M31Plmr   M11Plmr
## 549322   0   0   0       0       0       0 0.000e+00
## 522457   0   0   0       0       0       0 0.000e+00
## 951      0   0   0       0       0       0 2.305e-06
## 244423   0   0   0       0       0       0 0.000e+00
## 586076   0   0   0       0       0       0 0.000e+00
```


Here the dimension of the count data remains unchanged. This is in
fact a requirement for any `SummarizedExperiment` object.



### colData

`colData` contains data on the samples.


```r
colData(tse)
```

```
## DataFrame with 26 rows and 7 columns
##         X.SampleID   Primer Final_Barcode Barcode_truncated_plus_T
##           <factor> <factor>      <factor>                 <factor>
## CL3        CL3      ILBC_01        AACGCA                   TGCGTT
## CC1        CC1      ILBC_02        AACTCG                   CGAGTT
## SV1        SV1      ILBC_03        AACTGT                   ACAGTT
## M31Fcsw    M31Fcsw  ILBC_04        AAGAGA                   TCTCTT
## M11Fcsw    M11Fcsw  ILBC_05        AAGCTG                   CAGCTT
## ...            ...      ...           ...                      ...
## TS28         TS28   ILBC_25        ACCAGA                   TCTGGT
## TS29         TS29   ILBC_26        ACCAGC                   GCTGGT
## Even1        Even1  ILBC_27        ACCGCA                   TGCGGT
## Even2        Even2  ILBC_28        ACCTCG                   CGAGGT
## Even3        Even3  ILBC_29        ACCTGT                   ACAGGT
##         Barcode_full_length SampleType
##                    <factor>   <factor>
## CL3             CTAGCGTGCGT      Soil 
## CC1             CATCGACGAGT      Soil 
## SV1             GTACGCACAGT      Soil 
## M31Fcsw         TCGACATCTCT      Feces
## M11Fcsw         CGACTGCAGCT      Feces
## ...                     ...        ...
## TS28            GCATCGTCTGG      Feces
## TS29            CTAGTCGCTGG      Feces
## Even1           TGACTCTGCGG      Mock 
## Even2           TCTGATCGAGG      Mock 
## Even3           AGAGAGACAGG      Mock 
##                                        Description
##                                           <factor>
## CL3     Calhoun South Carolina Pine soil, pH 4.9  
## CC1     Cedar Creek Minnesota, grassland, pH 6.1  
## SV1     Sevilleta new Mexico, desert scrub, pH 8.3
## M31Fcsw M3, Day 1, fecal swab, whole body study   
## M11Fcsw M1, Day 1, fecal swab, whole body study   
## ...                                            ...
## TS28                                       Twin #1
## TS29                                       Twin #2
## Even1                                      Even1  
## Even2                                      Even2  
## Even3                                      Even3
```



### rowData

`rowData` contains data on the features of the analyzed samples. Of particular
interest for the microbiome field this is used to store taxonomic information.


```r
rowData(tse)
```

```
## DataFrame with 19216 rows and 7 columns
##            Kingdom        Phylum        Class        Order        Family
##        <character>   <character>  <character>  <character>   <character>
## 549322     Archaea Crenarchaeota Thermoprotei           NA            NA
## 522457     Archaea Crenarchaeota Thermoprotei           NA            NA
## 951        Archaea Crenarchaeota Thermoprotei Sulfolobales Sulfolobaceae
## 244423     Archaea Crenarchaeota        Sd-NA           NA            NA
## 586076     Archaea Crenarchaeota        Sd-NA           NA            NA
## ...            ...           ...          ...          ...           ...
## 278222    Bacteria           SR1           NA           NA            NA
## 463590    Bacteria           SR1           NA           NA            NA
## 535321    Bacteria           SR1           NA           NA            NA
## 200359    Bacteria           SR1           NA           NA            NA
## 271582    Bacteria           SR1           NA           NA            NA
##              Genus                Species
##        <character>            <character>
## 549322          NA                     NA
## 522457          NA                     NA
## 951     Sulfolobus Sulfolobusacidocalda..
## 244423          NA                     NA
## 586076          NA                     NA
## ...            ...                    ...
## 278222          NA                     NA
## 463590          NA                     NA
## 535321          NA                     NA
## 200359          NA                     NA
## 271582          NA                     NA
```

### rowTree  

Phylogenetic trees also play an important role for the microbiome field. The 
`TreeSummarizedExperiment` class is able to keep track of feature and node
relations via two functions, `rowTree` and `rowLinks`.

A tree can be accessed via `rowTree` as `phylo` object.       

```r
rowTree(tse)
```

```
## 
## Phylogenetic tree with 19216 tips and 19215 internal nodes.
## 
## Tip labels:
##   549322, 522457, 951, 244423, 586076, 246140, ...
## Node labels:
##   , 0.858.4, 1.000.154, 0.764.3, 0.995.2, 1.000.2, ...
## 
## Rooted; includes branch lengths.
```

The links to the individual features are available through `rowLinks`.


```r
rowLinks(tse)
```

```
## LinkDataFrame with 19216 rows and 5 columns
##           nodeLab   nodeNum nodeLab_alias    isLeaf   whichTree
##       <character> <integer>   <character> <logical> <character>
## 1          549322         1       alias_1      TRUE       phylo
## 2          522457         2       alias_2      TRUE       phylo
## 3             951         3       alias_3      TRUE       phylo
## 4          244423         4       alias_4      TRUE       phylo
## 5          586076         5       alias_5      TRUE       phylo
## ...           ...       ...           ...       ...         ...
## 19212      278222     19212   alias_19212      TRUE       phylo
## 19213      463590     19213   alias_19213      TRUE       phylo
## 19214      535321     19214   alias_19214      TRUE       phylo
## 19215      200359     19215   alias_19215      TRUE       phylo
## 19216      271582     19216   alias_19216      TRUE       phylo
```

Please note that there can be a 1:1 relationship between tree nodes and 
features, but this is not a must have. This means there can be features, which
are not linked to nodes, and nodes, which are not linked to features. To change
the links in an existing object, the `changeTree` function is available.



### Alternative experiments {#alt-exp}

_**Alternative experiments**_ differ from transformations as they can
contain complementary data, which is no longer tied to the same
dimensions as the assay data. However, the number of samples (columns)
must be the same.

This can come into play for instance when one has taxonomic abundance
profiles quantified with different measurement technologies, such as
phylogenetic microarrays, amplicon sequencing, or metagenomic
sequencing. Such alternative experiments that concern the same samples
can be stored as

1. Separate _assays_ assuming that the taxonomic information can be mapped 
between feature directly 1:1; or 
2. data in the _altExp_ slot of the `TreeSummarizedExperiment`, if the feature 
dimensions differ. Each element of the _altExp_ slot is a `SummarizedExperiment`
or an object from a derived class with independent feature data.


As an example, we show how to store taxonomic abundance tables
agglomerated at different taxonomic levels. However, the data could as
well originate from entirely different measurement sources as long as
the samples are matched.


```r
# Agglomerate the data to Phylym level
tse_phylum <- agglomerateByRank(tse, "Phylum")
# both have the same number of columns (samples)
dim(tse)
```

```
## [1] 19216    26
```

```r
dim(tse_phylum)
```

```
## [1] 67 26
```

```r
# Add the new table as an alternative experiment
altExp(tse, "Phylum") <- tse_phylum
altExpNames(tse)
```

```
## [1] "Phylum"
```

```r
# Pick a sample subset: this acts on both altExp and assay data
tse[,1:10]
```

```
## class: TreeSummarizedExperiment 
## dim: 19216 10 
## metadata(0):
## assays(2): counts relabundance
## rownames(19216): 549322 522457 ... 200359 271582
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(10): CL3 CC1 ... M31Tong M11Tong
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(1): Phylum
## rowLinks: a LinkDataFrame (19216 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```

```r
dim(altExp(tse[,1:10],"Phylum"))
```

```
## [1] 67 10
```

For more details of altExp have a look at the [Intro vignette](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) of the 
`SingleCellExperiment` package [@R-SingleCellExperiment].



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
`matrix`-like object including `SummarizedExperiment` objects, and the
number of samples can differ between the elements.



```r
#TODO: Find the right dataset to explain a non 1:1 sample relationship
```


For information have a look at the [intro vignette](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html) of the `MultiAssayExperiment` package.  

 
   Option   Rows (features)    Cols (samples)               Recommended  
---------   --------------    ---------------  ------------------------
   assays  	     match              match       Data transformations  
   altExp             free              match    Alternative experiments  
MultiAssay            free      free (mapping)    Multi-omic experiments    

Table: (\#tab:options) **Recommended options for storing multiple data tables in microbiome studies** The _assays_ are best suited for data transformations (one-to-one match between samples and columns across the assays). The _alternative experiments_ are particularly suitable for alternative versions of the data that are of same type but may have a different number of features (e.g. taxonomic groups); this is for instance the case with taxonomic abundance tables agglomerated at different levels (e.g. genus vs. phyla) or alternative profiling technologies (e.g. amplicon sequencing vs. shallow shotgun metagenomics). For alternative experiments one-to-one match between samples (cols) is required but the alternative experiment tables can have different numbers of features (rows). Finally, elements of the _MultiAssayExperiment_ provide the most flexible way to incorporate multi-omic data tables with flexible numbers of samples and features. We recommend these conventions as the basis for methods development and application in microbiome studies.




## Demonstration data {#example-data}

Open demonstration data for testing and benchmarking purposes is
available from multiple locations. This chapter introduces some
options. The other chapters of this book provide ample examples about
the use of the data.




### Package data {#package-data}

The `mia` R package contains example data sets that are direct
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


Check the documentation for this data set:





### ExperimentHub data

[ExperimentHub](https://bioconductor.org/packages/release/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html)
provides a variety of data resources, including the
[microbiomeDataSets](https://bioconductor.org/packages/release/data/experiment/html/microbiomeDataSets.html)
package [@Morgan2021; @microlahti2021].

A table of the available data sets is available through the
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
us retrieve a *[MultiAssayExperiment](https://bioconductor.org/packages/3.17/MultiAssayExperiment)* data set:


```r
# mae <- HintikkaXOData()
# Since HintikkaXOData is now added to mia, we can load it directly from there
# We suggest to check other datasets from microbiomeDataSets
data(HintikkaXOData)
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
is a large collection of curated human microbiome data sets, provided as
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

Result of amplicon sequencing is large number of files that include all the sequences
that were read from samples. Those sequences need to be matched with taxa. Additionally,
we need to know how many times each taxa were found from each sample. 

There are several algorithms to do that, and DADA2 is one of the most common. 
You can find DADA2 pipeline tutorial for example from 
[here](https://benjjneb.github.io/dada2/tutorial.html).
After DADA2 portion of the tutorial is the data is stored into _phyloseq_ object 
(Bonus: Handoff to phyloseq). To store the data to _TreeSummarizedExperiment_,
follow the example below. 

You can find full workflow script without further explanations and comments from 
[here](https://github.com/microbiome/OMA/blob/master/dada2_workflow.Rmd)



Load required packages.


```r
library(mia)
library(ggplot2)

if( !require("BiocManager") ){
    install.packages("BiocManager")
    library("BiocManager")
}

if( !require("Biostrings") ){
    BiocManager::install("Biostrings")
    library("Biostrings")
}
library(Biostrings)
```

Create arbitrary example sample metadata like it was done in tutorial. Usually, 
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

Convert data into right format and create _TreeSE_ object.


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



### Import from external files

Microbiome (taxonomic) profiling data is commonly distributed in
various file formats. You can import such external data files as a
(Tree)SummarizedExperiment object but the details depend on the file
format. Here, we provide examples for common formats. 


#### CSV import

**CSV data tables** can be imported with the standard R functions,
  then converted to the desired format. For detailed examples, you can
  check the [Bioconductor course
  material](https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html)
  by Martin Morgan. You can also check the [example
  files](https://github.com/microbiome/OMA/tree/master/data) and
  construct your own CSV files accordingly.

Recommendations for the CSV files are the following. File names are
arbitrary; we refer here to the same names than in the examples:

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
of the available demonstration data sets, and make sure that your data
components have the same format.

There are many different source files and many different ways to read
data in R. One can do data manipulation in R as well. Investigate the
entries as follows.



```r
# coldata rownames match assay colnames
all(rownames(samples) == colnames(counts)) # our data set
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
all(rownames(tax) == rownames(counts)) # our data set
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

This example shows how Biom files are imported into a
`TreeSummarizedExperiment` object.

The data is from following publication: 
Tengeler AC _et al._ (2020) [**Gut microbiota from persons with
attention-deficit/hyperactivity disorder affects the brain in
mice**](https://doi.org/10.1186/s40168-020-00816-x). 

The data set consists of 3 files:

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


## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.3.0 (2023-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] phyloseq_1.44.0                BiocManager_1.30.20           
 [3] ggplot2_3.4.2                  microbiomeDataSets_1.1.7      
 [5] mia_1.9.2                      MultiAssayExperiment_1.26.0   
 [7] TreeSummarizedExperiment_2.1.4 Biostrings_2.68.0             
 [9] XVector_0.40.0                 SingleCellExperiment_1.22.0   
[11] SummarizedExperiment_1.30.0    Biobase_2.60.0                
[13] GenomicRanges_1.52.0           GenomeInfoDb_1.36.0           
[15] IRanges_2.34.0                 S4Vectors_0.38.0              
[17] BiocGenerics_0.46.0            MatrixGenerics_1.12.0         
[19] matrixStats_0.63.0-9003        BiocStyle_2.28.0              
[21] rebook_1.9.0                  

loaded via a namespace (and not attached):
  [1] jsonlite_1.8.4                CodeDepends_0.6.5            
  [3] magrittr_2.0.3                ggbeeswarm_0.7.2             
  [5] rmarkdown_2.21                zlibbioc_1.46.0              
  [7] vctrs_0.6.2                   multtest_2.56.0              
  [9] memoise_2.0.1                 DelayedMatrixStats_1.22.0    
 [11] RCurl_1.98-1.12               htmltools_0.5.5              
 [13] AnnotationHub_3.8.0           curl_5.0.0                   
 [15] BiocNeighbors_1.18.0          Rhdf5lib_1.22.0              
 [17] rhdf5_2.44.0                  plyr_1.8.8                   
 [19] DECIPHER_2.28.0               cachem_1.0.7                 
 [21] igraph_1.4.2                  iterators_1.0.14             
 [23] mime_0.12                     lifecycle_1.0.3              
 [25] pkgconfig_2.0.3               rsvd_1.0.5                   
 [27] Matrix_1.5-4                  R6_2.5.1                     
 [29] fastmap_1.1.1                 GenomeInfoDbData_1.2.10      
 [31] shiny_1.7.4                   digest_0.6.31                
 [33] colorspace_2.1-0              AnnotationDbi_1.62.0         
 [35] scater_1.28.0                 irlba_2.3.5.1                
 [37] ExperimentHub_2.8.0           RSQLite_2.3.1                
 [39] vegan_2.6-4                   beachmat_2.16.0              
 [41] filelock_1.0.2                fansi_1.0.4                  
 [43] httr_1.4.5                    mgcv_1.8-42                  
 [45] compiler_4.3.0                withr_2.5.0                  
 [47] bit64_4.0.5                   BiocParallel_1.34.0          
 [49] viridis_0.6.2                 DBI_1.1.3                    
 [51] MASS_7.3-59                   rappdirs_0.3.3               
 [53] DelayedArray_0.25.0           biomformat_1.28.0            
 [55] permute_0.9-7                 tools_4.3.0                  
 [57] vipor_0.4.5                   beeswarm_0.4.0               
 [59] ape_5.7-1                     interactiveDisplayBase_1.38.0
 [61] httpuv_1.6.9                  glue_1.6.2                   
 [63] rhdf5filters_1.12.1           nlme_3.1-162                 
 [65] promises_1.2.0.1              grid_4.3.0                   
 [67] ade4_1.7-22                   cluster_2.1.4                
 [69] reshape2_1.4.4                generics_0.1.3               
 [71] gtable_0.3.3                  tidyr_1.3.0                  
 [73] data.table_1.14.8             BiocSingular_1.16.0          
 [75] ScaledMatrix_1.7.1            utf8_1.2.3                   
 [77] foreach_1.5.2                 ggrepel_0.9.3                
 [79] BiocVersion_3.17.1            pillar_1.9.0                 
 [81] stringr_1.5.0                 yulab.utils_0.0.6            
 [83] later_1.3.0                   splines_4.3.0                
 [85] dplyr_1.1.2                   BiocFileCache_2.8.0          
 [87] treeio_1.24.0                 lattice_0.21-8               
 [89] survival_3.5-5                bit_4.0.5                    
 [91] tidyselect_1.2.0              DirichletMultinomial_1.42.0  
 [93] scuttle_1.10.0                knitr_1.42                   
 [95] gridExtra_2.3                 bookdown_0.33                
 [97] xfun_0.39                     stringi_1.7.12               
 [99] lazyeval_0.2.2                yaml_2.3.7                   
[101] evaluate_0.20                 codetools_0.2-19             
[103] tibble_3.2.1                  graph_1.78.0                 
[105] cli_3.6.1                     xtable_1.8-4                 
[107] munsell_0.5.0                 Rcpp_1.0.10                  
[109] dir.expiry_1.8.0              dbplyr_2.3.2                 
[111] png_0.1-8                     XML_3.99-0.14                
[113] parallel_4.3.0                ellipsis_0.3.2               
[115] blob_1.2.4                    sparseMatrixStats_1.12.0     
[117] bitops_1.0-7                  decontam_1.20.0              
[119] viridisLite_0.4.1             tidytree_0.4.2               
[121] scales_1.2.1                  purrr_1.0.1                  
[123] crayon_1.5.2                  rlang_1.1.1                  
[125] KEGGREST_1.40.0              
```
</div>
