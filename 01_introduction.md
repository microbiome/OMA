# (PART) Introduction {-}

# Data Infrastructure {#data-introduction}

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

The
[`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
(`SE`) is a widely used class for analyzing data obtained by common sequencing
techniques. `SE` is common structure for several Bioconductor packages that are
used for analyzing RNAseq, ChIp-Seq data. `SE` class is also used in R packages
for analyzing microarrays, flow cytometry, proteomics, single-cell sequencing
data and many more. The single-cell analysis is facilitated by
[SingelCellExperiment
class](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)
(`SCE`), which allows the user to store results of dimensionality reduction and
alternative experiments. Alternative experiments (`altExps`) can be differently
processed data within the analysis workflows.

Recently,
[`TreeSummarizedExperiment`](https://www.bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html)
(`TSE`)
were developed to extend the `SE` and `SCE` class for incorporating hierarchical
information (including phylogenetic tree) and reference sequences.

The `mia` package implements tools using these classes for analysis of
microbiome sequencing data.

## Installation

Install the development version from GitHub using `remotes` R package.  


```r
# install remotes 
#install.packages("remotes")
BiocManager::install("FelixErnst/mia")
```

### Packages    

1. `mia`    : Microbiome analysis tools   
2. `miaViz` : Microbiome analysis specific visualization

**See also:**    

[`microbiome`](https://bioconductor.org/packages/devel/bioc/html/microbiome.html)


## Background

The widely used `phyloseq` package and class were around before the [`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)  
and the derived 
[`TreeSummarizedExperiment`](https://www.bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html) 
class. Many methods for taxonomic profiling data are readily for the  `phyloseq` class structure. 

In order to facilitate the transition, we provide here a short description how `phyloseq` and `*Experiment` classes relate to 
each other.

`assays`     : This slot is similar to `otu_table` in `phyloseq`. In a `SummarizedExperiment`
               object multiple assays, raw counts, transformed counts can be stored. See also 
               [`MultiAssayExperiment`](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html) for storing data from multiple `experiments` such as RNASeq, Proteomics, etc.       
`rowData`    : This slot is similar to `tax_table` in `phyloseq` to store taxonomic information.     
`colData`    : This slot is similar to `sample_data` in `phyloseq` to store information related to samples.    
`rowTree`    : This slot is similar to `phy_tree` in `phyloseq` to store phylogenetic tree.     

In this book, you will come across terms like `FeatureIDs` and `SampleIDs`.   
`FeatureIDs` : These are basically OTU/ASV ids which are row names in `assays` and `rowData`.    
`SampleIDs`  : As the name suggests, these are sample ids which are column names in `assays` and row names in `colData`.  

`FeatureIDs` and `SampleIDs` are used but the technical terms `rownames` and 
`colnames` are encouraged to be used, since they relate to actual objects we 
work with.

<img src="https://raw.githubusercontent.com/FelixErnst/TreeSummarizedExperiment
/2293440c6e70ae4d6e978b6fdf2c42fdea7fb36a/vignettes/tse2.png" width="100%"/>

**Figure sources:** 

**Original article**
-   Huang R _et al_. (2021) [TreeSummarizedExperiment: a S4 class 
for data with hierarchical structure](https://doi.org/10.12688/
f1000research.26669.2). F1000Research 9:1246.

**Reference Sequence slot extension**
- Lahti L _et al_. (2020) [Upgrading the R/Bioconductor ecosystem for microbiome 
research](https://doi.org/10.7490/
f1000research.1118447.1) F1000Research 9:1464 (slides).

## Loading experimental microbiome data


### Importing data from external files

Microbiome (taxonomic) profiling data is commonly distributed in
various file formats. You can import such external data files as a
(Tree)SummarizedExperiment object but the details depend on the file
format. Here, we provide examples for common formats.

**CSV data tables** can be imported with the standard R functions,
  then converted to the desired format. For detailed examples, you can
  check the [Bioconductor course
  material](https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html)
  by Martin Morgan. The following example reads abundance tables,
  taxonomic mapping tables, and sample metadata, assuming that the
  input data files are properly prepared with appropriate row and
  column names.


```r
counts <- read.csv(count_file)   # Abundance table (e.g. ASV data; to assay data)
tax <- read.csv(tax_file)        # Taxonomy table (to rowData)
samples <- read.csv(sample_file) # Sample data (to colData)
se <- SummarizedExperiment(assays = list(counts = counts),
                           colData = samples,
                           rowData = tax)
```

A specific import functions are provided for:

-   Biom files (see `help(mia::loadFromBiom)`)
-   QIIME2 files (see `help(mia::loadFromQIIME2)`)
-   Mothur files (see `help(mia::loadFromMothur)`)

#### Example for importing Biom files

This example shows how Biom files are imported into a TreeSummarizedExperiment object. 

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
## class: SummarizedExperiment 
## dim: 151 27 
## metadata(0):
## assays(1): counts
## rownames(151): 1726470 1726471 ... 17264756 17264757
## rowData names(6): taxonomy1 taxonomy2 ... taxonomy5 taxonomy6
## colnames(27): A110 A111 ... A38 A39
## colData names(0):
```

The `assays` slot includes a list of abundance tables. The imported abundance table is named as "counts".
Let us inspect only the first cols and rows.


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

These taxonomic rank names (column names) are not real rank names. Let’s replace them with real rank names.

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

Let us add a sample meta data file.


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

If the data is has already been imported in R in another format, it
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

# convert phyloseq to TSE
GlobalPatterns_TSE <- makeTreeSummarizedExperimentFromPhyloseq(GlobalPatterns_phyloseq) 
```


We can also convert `TreeSummarizedExperiment` objects into `phyloseq`
with respect to the shared components that are supported by both
formats (i.e. taxonomic abundance table, sample metadata, taxonomic
table, phylogenetic tree, sequence information). This is useful for
instance when additional methods are available for `phyloseq`.

TODO: conversion function from TSE to phyloseq





## Metadata

## Microbiome and tree data specific aspects


```r
library(mia)
data("GlobalPatterns")
se <- GlobalPatterns 
se
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

### Assays  

The `assays` slot contains the experimental data as count matrices. Multiple 
matrices can be stored the result of `assays` is actually a list of matrices.


```r
assays(se)
```

```
## List of length 1
## names(1): counts
```

Individual assays can be accessed via `assay`


```r
assay(se, "counts")[1:5,1:7]
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
calcualted and stored along the original count data using `relAbundanceCounts`.


```r
se <- relAbundanceCounts(se)
assays(se)
```

```
## List of length 2
## names(2): counts relabundance
```

Now there are two assays available in the `se` object, `counts` and 
`relabundance`.


```r
assay(se, "relabundance")[1:5,1:7]
```

```
##        CL3 CC1 SV1 M31Fcsw M11Fcsw M31Plmr   M11Plmr
## 549322   0   0   0       0       0       0 0.000e+00
## 522457   0   0   0       0       0       0 0.000e+00
## 951      0   0   0       0       0       0 2.305e-06
## 244423   0   0   0       0       0       0 0.000e+00
## 586076   0   0   0       0       0       0 0.000e+00
```

### colData

`colData` contains data on the samples.


```r
colData(se)
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

`rowData` contains data on the features of the samples analyzed. Of particular
interest for the microbiome field this is used to store taxonomic information.


```r
rowData(se)
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
rowTree(se)
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
rowLinks(se)
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

## Data conversion

Sometimes custom solutions are need for analyzing the data. `mia` contains a 
few functions to help in these situations.

### Tidy data

For several custom analysis and visualization packages such as those from the 
`tidyverse` the `SE` data can be converted to long data.frame format with 
`meltAssay`.    


```r
library(mia)
molten_se <- meltAssay(se,
                       add_row_data = TRUE,
                       add_col_data = TRUE,
                       abund_values = "relabundance")
molten_se
```

```
## # A tibble: 499,616 × 17
##    FeatureID SampleID counts Kingdom Phylum   Class   Order Family Genus Species
##    <fct>     <fct>     <dbl> <chr>   <chr>    <chr>   <chr> <chr>  <chr> <chr>  
##  1 549322    CL3           0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  2 549322    CC1           0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  3 549322    SV1           0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  4 549322    M31Fcsw       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  5 549322    M11Fcsw       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  6 549322    M31Plmr       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  7 549322    M11Plmr       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  8 549322    F21Plmr       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  9 549322    M31Tong       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
## 10 549322    M11Tong       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
## # … with 499,606 more rows, and 7 more variables: X.SampleID <fct>,
## #   Primer <fct>, Final_Barcode <fct>, Barcode_truncated_plus_T <fct>,
## #   Barcode_full_length <fct>, SampleType <fct>, Description <fct>
```

## Conclusion

Some wrapping up...

## Additional Reading

### Lecture slides

Introduction to microbiome data science [lecture slides](https://github.com/microbiome/course_2021_radboud/tree/main/slides).

### R programming resources

 * R programming basics: [Base R](https://www.rstudio.com/wp-content/uploads/2016/10/r-cheat-sheet-3.pdf)
 * Basics of R programming: [Base R](https://raw.githubusercontent.com/rstudio/cheatsheets/master/base-r.pdf)
 * [R cheat sheets](https://www.rstudio.com/resources/cheatsheets/)
 * R visualization with [ggplot2](https://www.rstudio.com/wp-content/uploads/2016/11/ggplot2-cheatsheet-2.1.pdf) 
 * [R graphics cookbook](http://www.cookbook-r.com/Graphs/)

Rmarkdown

* [Rmarkdown tips](https://rmarkdown.rstudio.com/)


RStudio

* [RStudio cheat sheet](https://www.rstudio.com/wp-content/uploads/2016/01/rstudio-IDE-cheatsheet.pdf) 

### Resources for TreeSummarizedExperiment

 * SingleCellExperiment
   + [Publication](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)
   + [Project page](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
 * SummarizedExperiment
   + [Publication](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html)
   + [Project page](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
 * TreeSummarizedExperiment
   + [Publication](https://f1000research.com/articles/9-1246)
   + [Project page](https://www.bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html)
   
### Resources for phyloseq

 * [List of R tools for microbiome analysis](https://microsud.github.io/Tools-Microbiome-Analysis/)
 * [phyloseq](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217)
 * [microbiome tutorial](http://microbiome.github.io/tutorials/)
 * [microbiomeutilities](https://microsud.github.io/microbiomeutilities/)
 * Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses ([Callahan et al. F1000, 2016](https://f1000research.com/articles/5-1492/v2)).




### Further reading


* [Data Analysis and Visualization in R for Ecologists](https://datacarpentry.org/R-ecology-lesson/) by Data Carpentry

* [Modern Statistics for Modern Biology. Holmes & Huber (2018)](http://web.stanford.edu/class/bios221/book/) for background in statistical analysis

* [Microbiome Data Science. Shetty & Lahti, 2019](https://openresearchlabs.github.io/publications/papers/2018-Shetty-Lahti-MDS.pdf)

## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.1.0 (2021-05-18)
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
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] mia_1.1.7                      TreeSummarizedExperiment_2.1.3
 [3] Biostrings_2.61.1              XVector_0.33.0                
 [5] SingleCellExperiment_1.15.1    SummarizedExperiment_1.23.1   
 [7] Biobase_2.53.0                 GenomicRanges_1.45.0          
 [9] GenomeInfoDb_1.29.3            IRanges_2.27.0                
[11] S4Vectors_0.31.0               BiocGenerics_0.39.1           
[13] MatrixGenerics_1.5.1           matrixStats_0.59.0            
[15] BiocStyle_2.21.3               rebook_1.3.0                  

loaded via a namespace (and not attached):
  [1] ggbeeswarm_0.6.0            colorspace_2.0-2           
  [3] ellipsis_0.3.2              scuttle_1.3.0              
  [5] BiocNeighbors_1.11.0        bit64_4.0.5                
  [7] fansi_0.5.0                 decontam_1.13.0            
  [9] splines_4.1.0               codetools_0.2-18           
 [11] sparseMatrixStats_1.5.0     cachem_1.0.5               
 [13] knitr_1.33                  scater_1.21.2              
 [15] ade4_1.7-17                 phyloseq_1.37.0            
 [17] jsonlite_1.7.2              cluster_2.1.2              
 [19] graph_1.71.2                BiocManager_1.30.16        
 [21] compiler_4.1.0              assertthat_0.2.1           
 [23] Matrix_1.3-4                fastmap_1.1.0              
 [25] lazyeval_0.2.2              cli_3.0.1                  
 [27] BiocSingular_1.9.1          htmltools_0.5.1.1          
 [29] tools_4.1.0                 igraph_1.2.6               
 [31] rsvd_1.0.5                  gtable_0.3.0               
 [33] glue_1.4.2                  GenomeInfoDbData_1.2.6     
 [35] reshape2_1.4.4              dplyr_1.0.7                
 [37] Rcpp_1.0.7                  jquerylib_0.1.4            
 [39] rhdf5filters_1.5.0          vctrs_0.3.8                
 [41] multtest_2.49.0             ape_5.5                    
 [43] nlme_3.1-152                DECIPHER_2.21.0            
 [45] iterators_1.0.13            DelayedMatrixStats_1.15.0  
 [47] xfun_0.24                   stringr_1.4.0              
 [49] beachmat_2.9.0              lifecycle_1.0.0            
 [51] irlba_2.3.3                 XML_3.99-0.6               
 [53] zlibbioc_1.39.0             MASS_7.3-54                
 [55] scales_1.1.1                biomformat_1.21.0          
 [57] parallel_4.1.0              rhdf5_2.37.0               
 [59] yaml_2.2.1                  memoise_2.0.0              
 [61] gridExtra_2.3               ggplot2_3.3.5              
 [63] sass_0.4.0                  stringi_1.7.3              
 [65] RSQLite_2.2.7               foreach_1.5.1              
 [67] ScaledMatrix_1.1.0          tidytree_0.3.4             
 [69] permute_0.9-5               filelock_1.0.2             
 [71] BiocParallel_1.27.2         rlang_0.4.11               
 [73] pkgconfig_2.0.3             bitops_1.0-7               
 [75] evaluate_0.14               lattice_0.20-44            
 [77] Rhdf5lib_1.15.2             purrr_0.3.4                
 [79] treeio_1.17.2               CodeDepends_0.6.5          
 [81] bit_4.0.4                   tidyselect_1.1.1           
 [83] plyr_1.8.6                  magrittr_2.0.1             
 [85] bookdown_0.22               R6_2.5.0                   
 [87] generics_0.1.0              DelayedArray_0.19.1        
 [89] DBI_1.1.1                   mgcv_1.8-36                
 [91] pillar_1.6.1                survival_3.2-11            
 [93] RCurl_1.98-1.3              tibble_3.1.3               
 [95] dir.expiry_1.1.0            crayon_1.4.1               
 [97] utf8_1.2.2                  rmarkdown_2.9              
 [99] viridis_0.6.1               grid_4.1.0                 
[101] data.table_1.14.0           blob_1.2.2                 
[103] vegan_2.5-7                 digest_0.6.27              
[105] tidyr_1.1.3                 munsell_0.5.0              
[107] DirichletMultinomial_1.35.0 beeswarm_0.4.0             
[109] viridisLite_0.4.0           vipor_0.4.5                
[111] bslib_0.2.5.1              
```
</div>
