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

[`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
(`SE`) is a generic and highly optimized container for complex data
structures. It has become a common choice for analysing various types
of biomedical profiling data, such as RNAseq, ChIp-Seq, microarrays,
flow cytometry, proteomics, and single-cell
sequencing. [`TreeSummarizedExperiment`](https://www.bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html)
(`TreeSE`) was developed as an extension to incorporate hierarchical
information (such as phylogenetic trees and sample hierarchies) and
reference sequences.

In microbiome data science, these containers can be used to link
taxonomic abundance tables with rich side information on the features
and samples. Taxonomic abundance data can be obtained by 16S rRNA
amplicon or metagenomic sequencing, phylogenetic microarrays, or by
other means.



## Package installation

Several R packages are now available providing methods for the
analysis and manipulation of data sets in such containers. One of
these is `mia`.

Install the Biocondcuctor development version with:


```r
BiocManager::install("microbiome/mia", version="devel")
```


## Background

The widely used `phyloseq` package and class were around before the
[`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
and the derived
[`TreeSummarizedExperiment`](https://www.bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html)
class. Many methods for taxonomic profiling data are readily available
for the `phyloseq` class structure.

In order to facilitate the transition, we provide here a short
description how `phyloseq` and `*Experiment` classes relate to each
other.

`assays` : This slot is similar to `otu_table` in `phyloseq`. In a
               `SummarizedExperiment` object multiple assays, raw
               counts, transformed counts can be stored. See also
               [`MultiAssayExperiment`](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html)
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
counts  <- read.csv(count_file)   # Abundance table (e.g. ASV data; to assay data)
tax     <- read.csv(tax_file)     # Taxonomy table (to rowData)
samples <- read.csv(sample_file)  # Sample data (to colData)
se <- SummarizedExperiment(assays = list(counts = counts),
                           colData = samples,
                           rowData = tax)
```

Specific import functions are provided for:

-   Biom files (see `help(mia::loadFromBiom)`)
-   QIIME2 files (see `help(mia::loadFromQIIME2)`)
-   Mothur files (see `help(mia::loadFromMothur)`)


#### Example for importing Biom files

This example shows how Biom files are imported into a `TreeSummarizedExperiment` object. 

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

## Metadata

## Microbiome and tree data specific aspects


```r
library(mia)
data("GlobalPatterns", package = "mia")
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

`rowData` contains data on the features of the analyzed samples. Of particular
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

## Data manipulation

Sometimes custom solutions are needed for analyzing the data. `mia` contains a 
few functions to help in these situations.

### Tidy data

For several custom analysis and visualization packages, such as those from the 
`tidyverse`, the `SE` data can be converted to long data.frame format with 
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
##    FeatureID SampleID relabundance Kingdom Phylum    Class    Order Family Genus
##    <fct>     <fct>           <dbl> <chr>   <chr>     <chr>    <chr> <chr>  <chr>
##  1 549322    CL3                 0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  2 549322    CC1                 0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  3 549322    SV1                 0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  4 549322    M31Fcsw             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  5 549322    M11Fcsw             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  6 549322    M31Plmr             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  7 549322    M11Plmr             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  8 549322    F21Plmr             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  9 549322    M31Tong             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
## 10 549322    M11Tong             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
## # … with 499,606 more rows, and 8 more variables: Species <chr>,
## #   X.SampleID <fct>, Primer <fct>, Final_Barcode <fct>,
## #   Barcode_truncated_plus_T <fct>, Barcode_full_length <fct>,
## #   SampleType <fct>, Description <fct>
```

### Subsetting

**Subsetting** data helps to draw the focus of analysis on particular
  sets of samples and / or features. When dealing with large data
  sets, the subset of interest can be extracted and investigated
  separately. This might improve performance and reduce the
  computational load.

Load:

* mia
* dplyr
* knitr
* data `GlobalPatterns`



Let us store `GlobalPatterns` into `se` and check its original number of features (rows) and samples (columns). **Note**: when subsetting by sample, expect the number of columns to decrease; when subsetting by feature, expect the number of rows to decrease.


```r
# store data into se and check dimensions
data("GlobalPatterns", package="mia")
se <- GlobalPatterns
# show dimensions
dim(se) ## [1] n°features    n°samples
```

```
## [1] 19216    26
```

#### Subset by sample (column-wise)

For the sake of demonstration, here we will extract a subset containing only the samples of human origin (feces, skin or tongue), stored as `SampleType` within `colData(se)` and also in `se`.

First, we would like to see all the possible values that `SampleType` can take on and how frequent those are: 


```r
# inspect possible values for SampleType
unique(se$SampleType)
```

```
## [1] Soil               Feces              Skin               Tongue            
## [5] Freshwater         Freshwater (creek) Ocean              Sediment (estuary)
## [9] Mock              
## 9 Levels: Feces Freshwater Freshwater (creek) Mock ... Tongue
```

```r
# show recurrence for each value
se$SampleType %>% table()
```
<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> . </th>
   <th style="text-align:right;"> Freq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Feces </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Freshwater </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Freshwater (creek) </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mock </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ocean </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sediment (estuary) </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Skin </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Soil </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Tongue </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
</tbody>
</table></div>

**Note**: after subsetting, expect the number of columns to equal the
  sum of the recurrences of the samples that you are interested
  in. For instance, `ncols = Feces + Skin + Tongue = 4 + 3 + 2 = 9`.

Next, we _logical index_ across the columns of `se` (make sure to
leave the first index empty to select all rows) and filter for the
samples of human origin. For this, we use the information on the
samples from the meta data `colData(se)`.


```r
# subset by sample
se_subset_by_sample <- se[ , se$SampleType %in% c("Feces", "Skin", "Tongue")]

# show dimensions
dim(se_subset_by_sample)
```

```
## [1] 19216     9
```

As a sanity check, the new object `se_subset_by_sample` should have
the original number of features (rows) and a number of samples
(columns) equal to the sum of the samples of interest (in this case
9).

Several characteristics can be used to subset by sample:

* origin
* sampling time
* sequencing method
* DNA / RNA barcode
* cohort

#### Subset by feature (row-wise)

Similarly, here we will extract a subset containing only the features
that belong to the Phyla "Actinobacteria" and "Chlamydiae", stored as
`Phylum` within `rowData(se)`. However, subsetting by feature implies
a few more obstacles, such as the presence of NA elements and the
possible need for agglomeration.

As previously, we would first like to see all the possible values that
`Phylum` can take on and how frequent those are:
  

```r
# inspect possible values for Phylum
unique(rowData(se)$Phylum)
```

```
##  [1] "Crenarchaeota"    "Euryarchaeota"    "Actinobacteria"   "Spirochaetes"    
##  [5] "MVP-15"           "Proteobacteria"   "SBR1093"          "Fusobacteria"    
##  [9] "Tenericutes"      "ZB3"              "Cyanobacteria"    "GOUTA4"          
## [13] "TG3"              "Chlorobi"         "Bacteroidetes"    "Caldithrix"      
## [17] "KSB1"             "SAR406"           "LCP-89"           "Thermi"          
## [21] "Gemmatimonadetes" "Fibrobacteres"    "GN06"             "AC1"             
## [25] "TM6"              "OP8"              "Elusimicrobia"    "NC10"            
## [29] "SPAM"             NA                 "Acidobacteria"    "CCM11b"          
## [33] "Nitrospirae"      "NKB19"            "BRC1"             "Hyd24-12"        
## [37] "WS3"              "PAUC34f"          "GN04"             "GN12"            
## [41] "Verrucomicrobia"  "Lentisphaerae"    "LD1"              "Chlamydiae"      
## [45] "OP3"              "Planctomycetes"   "Firmicutes"       "OP9"             
## [49] "WPS-2"            "Armatimonadetes"  "SC3"              "TM7"             
## [53] "GN02"             "SM2F11"           "ABY1_OD1"         "ZB2"             
## [57] "OP11"             "Chloroflexi"      "SC4"              "WS1"             
## [61] "GAL15"            "AD3"              "WS2"              "Caldiserica"     
## [65] "Thermotogae"      "Synergistetes"    "SR1"
```

```r
# show recurrence for each value
rowData(se)$Phylum %>% table()
```
<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> . </th>
   <th style="text-align:right;"> Freq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ABY1_OD1 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AC1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acidobacteria </td>
   <td style="text-align:right;"> 1021 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Actinobacteria </td>
   <td style="text-align:right;"> 1631 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AD3 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Armatimonadetes </td>
   <td style="text-align:right;"> 61 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bacteroidetes </td>
   <td style="text-align:right;"> 2382 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BRC1 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caldiserica </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caldithrix </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCM11b </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chlamydiae </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chlorobi </td>
   <td style="text-align:right;"> 64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chloroflexi </td>
   <td style="text-align:right;"> 437 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Crenarchaeota </td>
   <td style="text-align:right;"> 106 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cyanobacteria </td>
   <td style="text-align:right;"> 393 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Elusimicrobia </td>
   <td style="text-align:right;"> 31 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Euryarchaeota </td>
   <td style="text-align:right;"> 102 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fibrobacteres </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Firmicutes </td>
   <td style="text-align:right;"> 4356 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fusobacteria </td>
   <td style="text-align:right;"> 37 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAL15 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gemmatimonadetes </td>
   <td style="text-align:right;"> 191 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GN02 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GN04 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GN06 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GN12 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GOUTA4 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyd24-12 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KSB1 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCP-89 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LD1 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lentisphaerae </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MVP-15 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NC10 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nitrospirae </td>
   <td style="text-align:right;"> 74 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NKB19 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OP11 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OP3 </td>
   <td style="text-align:right;"> 30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OP8 </td>
   <td style="text-align:right;"> 27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OP9 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PAUC34f </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Planctomycetes </td>
   <td style="text-align:right;"> 638 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Proteobacteria </td>
   <td style="text-align:right;"> 6416 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAR406 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SBR1093 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SC3 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SC4 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SM2F11 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SPAM </td>
   <td style="text-align:right;"> 22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Spirochaetes </td>
   <td style="text-align:right;"> 124 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SR1 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synergistetes </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Tenericutes </td>
   <td style="text-align:right;"> 143 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TG3 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermi </td>
   <td style="text-align:right;"> 46 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermotogae </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TM6 </td>
   <td style="text-align:right;"> 27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TM7 </td>
   <td style="text-align:right;"> 32 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Verrucomicrobia </td>
   <td style="text-align:right;"> 470 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WPS-2 </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WS1 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WS2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WS3 </td>
   <td style="text-align:right;"> 70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZB2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZB3 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
</tbody>
</table></div>

**Note**: after subsetting, expect the number of columns to equal the
  sum of the recurrences of the feature(s) that you are interested
  in. For instance, `nrows = Actinobacteria + Chlamydiae = 1631 + 21 =
  1652`.

Depending on your research question, you might need to or need not
agglomerate the data in the first place: if you want to find the
abundance of each and every feature that belongs to Actinobacteria and
Chlamydiae, agglomeration is not needed; if you want to find the total
abundance of all the features that belong to Actinobacteria or
Chlamydiae, agglomeration is recommended.

##### Non-agglomerated data

Next, we _logical index_ across the rows of `se` (make sure to leave
the second index empty to select all columns) and filter for the
features that fall in either Actinobacteria or Chlamydiae. For this,
we use the information on the samples from the meta data
`rowData(se)`.

The first term with the `%in%` operator are includes all the features
of interest, whereas the second term after the AND operator `&`
filters out all the features that present a NA in place of Phylum.


```r
# subset by feature
se_subset_by_feature <- se[rowData(se)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(se)$Phylum), ]

# show dimensions
dim(se_subset_by_feature)
```

```
## [1] 1652   26
```

As a sanity check, the new object `se_subset_by_feature` should have the original number of samples (columns) and a number of features (rows) equal to the sum of the features of interest (in this case 1652).

##### Agglomerated data

When total abundances of certain Phyla are of relevance, the data is initially agglomerated by Phylum. Then, similar steps as in the case of not agglomerated data are followed.


```r
# agglomerate by Phylum
se_phylum <- se %>% agglomerateByRank(rank = "Phylum")

# subset by feature and get rid of NAs
se_phylum_subset_by_feature <- se_phylum[rowData(se_phylum)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(se_phylum)$Phylum), ]

# show dimensions
dim(se_phylum_subset_by_feature)
```

```
## [1]  2 26
```

**Note**: as data was agglomerated, the number of rows equal the
  number of Phyla used to index (in this case, just 2)

Alternatively:


```r
# store features of interest into phyla
phyla <- c("Phylum:Actinobacteria", "Phylum:Chlamydiae")
# subset by feature
se_phylum_subset_by_feature <- se_phylum[phyla, ]
# show dimensions
dim(se_subset_by_feature)
```

```
## [1] 1652   26
```

The code above returns the not agglomerated version of the data.

Fewer characteristics can be used to subset by feature:

* Taxonomic rank
* Meta-taxonomic group

For subsetting by Kingdom, agglomeration does not apply, whereas for
the other ranks it can be applied if necessary.

#### Subset by sample and feature

Finally, we can subset data by sample and feature at once. The
resulting subset contains all the samples of human origin and all the
features of Phyla "Actinobacteria" or "Chlamydiae".


```r
# subset by sample and feature and get rid of NAs
se_subset_by_sample_feature <- se[rowData(se)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(se)$Phylum), se$SampleType %in% c("Feces", "Skin", "Tongue")]

# show dimensions
dim(se_subset_by_sample_feature)
```

```
## [1] 1652    9
```

**Note**: the dimensions of `se_subset_by_sample_feature` agree with
  those of the previous subsets (9 columns filtered by sample and 1652
  rows filtered by feature).

If a study was to consider and quantify the presence of Actinobacteria
as well as Chlamydiae in different sites of the human body,
`se_subset_by_sample_feature` might be a suitable subset to start
with.


## Additional Reading


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
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] knitr_1.33                     dplyr_1.0.7                   
 [3] phyloseq_1.37.0                mia_1.1.14                    
 [5] TreeSummarizedExperiment_2.1.4 Biostrings_2.61.2             
 [7] XVector_0.33.0                 SingleCellExperiment_1.15.2   
 [9] SummarizedExperiment_1.23.4    Biobase_2.53.0                
[11] GenomicRanges_1.45.0           GenomeInfoDb_1.29.8           
[13] IRanges_2.27.2                 S4Vectors_0.31.3              
[15] BiocGenerics_0.39.2            MatrixGenerics_1.5.4          
[17] matrixStats_0.60.1-9001        BiocStyle_2.21.3              
[19] rebook_1.3.1                  

loaded via a namespace (and not attached):
  [1] systemfonts_1.0.2           plyr_1.8.6                 
  [3] igraph_1.2.6                lazyeval_0.2.2             
  [5] splines_4.1.1               BiocParallel_1.27.6        
  [7] ggplot2_3.3.5               scater_1.21.3              
  [9] digest_0.6.27               foreach_1.5.1              
 [11] htmltools_0.5.2             viridis_0.6.1              
 [13] fansi_0.5.0                 magrittr_2.0.1             
 [15] memoise_2.0.0               ScaledMatrix_1.1.0         
 [17] cluster_2.1.2               DECIPHER_2.21.0            
 [19] svglite_2.0.0               colorspace_2.0-2           
 [21] blob_1.2.2                  rvest_1.0.1                
 [23] ggrepel_0.9.1               xfun_0.25                  
 [25] crayon_1.4.1                RCurl_1.98-1.4             
 [27] jsonlite_1.7.2              graph_1.71.2               
 [29] survival_3.2-13             iterators_1.0.13           
 [31] ape_5.5                     glue_1.4.2                 
 [33] kableExtra_1.3.4            gtable_0.3.0               
 [35] zlibbioc_1.39.0             webshot_0.5.2              
 [37] DelayedArray_0.19.2         BiocSingular_1.9.1         
 [39] Rhdf5lib_1.15.2             scales_1.1.1               
 [41] DBI_1.1.1                   Rcpp_1.0.7                 
 [43] viridisLite_0.4.0           decontam_1.13.0            
 [45] tidytree_0.3.5              bit_4.0.4                  
 [47] rsvd_1.0.5                  httr_1.4.2                 
 [49] dir.expiry_1.1.0            ellipsis_0.3.2             
 [51] pkgconfig_2.0.3             XML_3.99-0.7               
 [53] scuttle_1.3.1               CodeDepends_0.6.5          
 [55] sass_0.4.0                  utf8_1.2.2                 
 [57] tidyselect_1.1.1            rlang_0.4.11               
 [59] reshape2_1.4.4              munsell_0.5.0              
 [61] tools_4.1.1                 cachem_1.0.6               
 [63] cli_3.0.1                   DirichletMultinomial_1.35.0
 [65] generics_0.1.0              RSQLite_2.2.8              
 [67] ade4_1.7-17                 evaluate_0.14              
 [69] biomformat_1.21.0           stringr_1.4.0              
 [71] fastmap_1.1.0               yaml_2.2.1                 
 [73] bit64_4.0.5                 purrr_0.3.4                
 [75] nlme_3.1-153                sparseMatrixStats_1.5.3    
 [77] xml2_1.3.2                  compiler_4.1.1             
 [79] rstudioapi_0.13             beeswarm_0.4.0             
 [81] filelock_1.0.2              treeio_1.17.2              
 [83] tibble_3.1.4                bslib_0.3.0                
 [85] stringi_1.7.4               highr_0.9                  
 [87] lattice_0.20-44             Matrix_1.3-4               
 [89] vegan_2.5-7                 permute_0.9-5              
 [91] multtest_2.49.0             vctrs_0.3.8                
 [93] pillar_1.6.2                lifecycle_1.0.0            
 [95] rhdf5filters_1.5.0          BiocManager_1.30.16        
 [97] jquerylib_0.1.4             BiocNeighbors_1.11.0       
 [99] data.table_1.14.0           bitops_1.0-7               
[101] irlba_2.3.3                 R6_2.5.1                   
[103] bookdown_0.24               gridExtra_2.3              
[105] vipor_0.4.5                 codetools_0.2-18           
[107] MASS_7.3-54                 assertthat_0.2.1           
[109] rhdf5_2.37.3                GenomeInfoDbData_1.2.6     
[111] mgcv_1.8-36                 parallel_4.1.1             
[113] grid_4.1.1                  beachmat_2.9.1             
[115] tidyr_1.1.3                 rmarkdown_2.10             
[117] DelayedMatrixStats_1.15.4   ggbeeswarm_0.6.0           
```
</div>
