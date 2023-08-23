# Exercises {#exercises}

Here you can find assignments on different topics. 

**Tips for exercises:**

   - Add comments that explain what each line or lines of code do. This helps you and others understand your code and find bugs. Furthermore, it is easier for you to reuse the code, and it promotes transparency.
   - Interpret results by comparing them to literature. List main findings, so that results can easily be understood by others without advanced knowledge on data analytics.
   - Avoid hard-coding. Use variables which get values in the beginning of the pipeline. That way it is easier for you to modify parameters and reuse the code.

## Workflows

### Reproducible reporting with Quarto

The following batch of exercises walks you through typical use cases of Quarto
in RStudio. Before heading to the exercises, it is recommended to read the
[Quarto guidelines for RStudio](https://quarto.org/docs/tools/rstudio.html)

#### New document

This exercise gets you started with creating a Quarto document and adding
text to it with typing conventions borrowed from the
[markdown syntax](https://quarto.org/docs/authoring/markdown-basics.html).
Feel free to render the document with the **Render** button after each step to
see the changes in the final report.

1. Open RStudio and create a new Quarto file named `My first Quarto`.
2. Add the subtitle `My first section` and write some text of your choice
   underneath. You can choose the level of headings by the number of preceding
   hashes (`#`).
3. Add a subsection named `List of items` and list three items underneath, both
   ordered and unordered. You can initialize items with numbers (`1.`, `2.`,
   `3.`, ...) or stars (`*`) for the ordered and unordered case, respectively.
4. Add another subsection named `Link to web` and add a clickable link to the
   [OMA book](https://microbiome.github.io/OMA/), using the pattern `[text](url)`.
5. Render the document and check its appearance

Nice start! You are now able to create a Quarto document, understand its syntax
and can render it into a reproducible report. If you got stuck, you can look up
the docs on [creating](https://quarto.org/docs/tools/rstudio.html#creating-documents)
and [rendering](https://quarto.org/docs/tools/rstudio.html#render-and-preview)
Quarto documents.

#### Code chunks

While customizable text is nothing new by itself, the advantage of Quarto (and
previously Rmakdown) is to combine text with code in R or other programming
languages, so that both the analytical pipeline and verbal description can be
put in one place. In this exercise, we learn how to write and run code in Quarto.

1. Open RStudio and create a new Quarto file.
2. Initialize a code chunk by pressing `alt` + `cmd` + `i` and define the variables
   `A <- "my name"` and `B <- 0` in it.
3. Write the text `Below is my first code chunk` just above the code chunk.
4. Initialize one more code chunk and add 100 to the variable `B` in it. 
5. Write the text `Below I change variable B` just above the second chunk.
6. **Extra**: Write the following line of text: `my name is A and I am B years old`,
   where A and B are variables defined in the code chunks upstream and change if
   those variables are modified. Be aware that inline code can be added as
   `> r my_inline_command` (without `>`).
   
Good job. You are now able to combine text and code in a Quarto document. If you
got stuck, you can refer to the Quarto docs on
[using code chunks](https://quarto.org/docs/visual-editor/technical.html#code-chunks).

#### Knitr options

Code chunks can be greatly customized in terms of visibility and execution,
output size and location and much more. This is possible with the knitr chunk
options, which usually appear at the beginning of the target chunk with the
syntax `#| option-name: value`, also described
[here](https://quarto.org/docs/computations/r.html#chunk-options).
In this exercise, we explore part of the knitr potential.

1. Open RStudio and create a new Quarto file.
2. Initialize three code chunks and label them as `setup`, `fig-box` and
   `tbl-coldata`, respectively. Remember that the name of a chunk can be
   specified with the `label` option.
3. Write the following code in the corresponding chunk and render the document.


```r
# setup
library(mia)
data("GlobalPatterns", package = "mia")
tse <- GlobalPatterns

# this line sets some options for all the chunks (global chunk options)
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
```


```r
# fig-box
boxplot(colSums(assay(tse)) ~ tse$SampleType)
```


```r
# tbl-coldata
knitr::kable(head(colData(tse)))
```

4. Set `include: false` in the `setup` chunk, `fig-width: 10` in the `fig-box`
   chunk and `echo: false` in the `tbl-coldata` chunk. Render the document again
   and find the differences from before.
5. Add the options `fig-cap` and `tab-cap` to the `fig-box` and `tbl-coldata`
   chunks, respectively. They require some text input, which makes for the
   caption of the figures or tables.
6. **Extra**: Create a cross-reference to `fig-box` and `tbl-coldata` in the
   text above the respective code chunk. You can do that with the syntax
   `@chunk-name`.
7. **Extra**: Define a custom folder for storing figures with `fig-path`. Insert
   it in `knitr::opts_chunk$set`, so that it applies globally to all the figures
   generated in the document.
   
Congratulations! You are now familiar with the great potential and flexibility
of knitr chunk options. An exhaustive list of available options can be found
in the [knitr documentation](https://yihui.org/knitr/options/).

#### YAML instructions

The box at the beginning of every Quarto document contains yaml options that let
you define the metadata of the document. They will affect the appearance of the
document when it is rendered. By default, the box includes yaml options for the
title, format and editor to be used, but much more information on layout, code
execution and figures can be specified. A comprehensive list of yaml options is
available [here](https://quarto.org/docs/reference/formats/html.html). In this
exercise, we will get a tiny taste of such functionality.

1. Open RStudio and create a new Quarto file.
2. In the yaml box at the beginning of the document, change the title from
   `Untitled` to `My first Quarto`.
4. In the same box, add the two lines `author` and `date` followed by your name
   and today's date, respectively.
5. Render the document and check its appearance.
6. **Extra**: Set `toc: true` to produce a table of contents. This line should
   follow `format` and `html` at the second level of indentation.
   
Well done! Now you are able to specify yaml options and understand how they
affect your Quarto document. If you got stuck, you can check
[this section](https://quarto.org/docs/tools/rstudio.html#yaml-intelligence) of
the Quarto documentation.

#### Quarto parameters

An advanced feature of Quarto consists of execution parameters, which are
externally pre-defined variables that are also accessible in the Quarto document.
They can be specified in the yaml box as `params`. Here we learn how to use them.

1. Open RStudio and create a new Quarto file.
2. In the yaml box at the beginning of the document, add a line named `params`
   followed by an indented line with `gamma: 10`
3. Initialize a code chunk and type `str(params$gamma)` in it.
4. Render the document and check what happened.
5. Define one more parameter `beta: 3` and multiply `gamma` by `beta` in a code
   chunk below.
6. Render the document again and check what happened.

Well done! You can now use an advanced feature of Quarto such as parameters.
If you got stuck, [here](https://quarto.org/docs/computations/parameters.html#knitr)
you can find more information about parameter definition and usage.


## Data containers: TreeSE

TreeSE containers represent the working unit of the mia package. In the
following exercises we learn how to construct, explore and work with them.
A few demo datasets can be imported with mia and can be accessed as explained
in chapter \@ref(example-data).

### Constructing a data object {#construct-TreeSE}

Here we cover how to construct a TreeSE from CSV files, using the components
of OKeefeDSData from the microbiomeDataSets package as an example dataset.

1. Fetch or download the files in [this directory](https://github.com/microbiome/data/tree/main/OKeefeDSData).
2. Read in the csv files with `read.csv` and store them into the variables
   `assays`, `rowdata` and `coldata`, respectively.
3. Create a TreeSE from the individual components with
   `TreeSummarizedExperiment`. Note that the function requires three arguments:
   assays, rowData and colData, to which you can give the appropriate item.
4. Check that importing is done correctly. E.g., choose random samples and
   features, and check that their values equal between raw files and TreeSE.

Usefuls functions: DataFrame, TreeSummarizedExperiment, matrix, rownames, colnames, SimpleList

### Importing data

Raw data of different types can be imported as a TreeSE with a number of
functions explained in chapter \@ref(#import-from-file). You can also check the
[function reference in the mia package](https://microbiome.github.io/mia/reference/index.html).

1. Get familiar with the
   [microbiome data repository](https://github.com/microbiome/data)
   and read the instructions in its
   [README](https://github.com/microbiome/data#training-microbiome-datasets)
   to import and construct datasets from there.
2. Import data from another format (functions: loadFromMetaphlan | loadFromMothur | loadFromQIIME2 | makeTreeSummarizedExperimentFromBiom | makeTreeSummarizedExperimentFromDADA2 ...)
3. Try out conversions between TreeSE and phyloseq data containers (makeTreeSummarizedExperimentFromPhyloseq; makephyloseqFromTreeSummarizedExperiment)

### Preliminary exploration

1. Import the mia package, load peerj13075 with `data` and store it into a
   variable named `tse`.
2. Get a summary about the TreeSE with `summary`. What is the mean count across
   samples? How many features recur only once (singletons)?
3. Check the dimensions of the TreeSE with `dim` or alternatively with `nrow`
   and `ncol`. How many samples and features are present?
4. List sample and features names with `rownames` and `colnames`.
5. Check what information about samples and features is contained by the
   `colData` and `rowData` of the TreeSE with `names`.
5. **Extra**: Calculate the number of unique taxa for each taxonomic rank. You
   can use `apply` to count unique elements for each column of rowData.

### Assay retrieval

1. Import the mia package, load peerj13075 with `data` and store it into a
   variable named `tse`.
2. List the names of all available assays with `assayNames`.
3. Fetch the list of assays with `assays`.
4. Retrieve the first assay of the TreeSE with `assay`, where the second
   argument can be either the name or the index of the desired assay.
   
Well done! You can now locate and retrieve individual assays of a TreeSE. If you
got stuck, you can refer to chapter \@ref(assay-slot) of this book.


### Sample information 

1. Import the mia package, load peerj13075 with `data` and store it into a
   variable named `tse`.
2. Check the names of the samples with `colnames`.
3. List the information on samples available in colData with `names`.
4. Visualize the colData with `View` and briefly look at the information stored
   in the different columns.
3. Get the abundances of all features for a specific sample, such as `ID34`. You
   can access sample-specific abundances with `getAbundaceSample`, specifying the
   desired `sample_id` and `assay.type`.


### Feature information

1. Import the mia package, load peerj13075 with `data` and store it into a
   variable named `tse`.
2. Check the names of the features with `rownames`.
3. List the information on features available in rowData with `names`.
4. Visualize the rowData with `View` and briefly look at the information stored
   in the different columns.
5. Get the abundances for a specific feature, such as `OTU1810`, in all the
   samples. You can access feature-specific abundances with `getAbundaceFeature`,
   specifying the desired `feature_id` and `assay.type`.
6. **Extra**: Create a taxonomy tree based on the taxonomy mappings with
   `addTaxonomyTree` and display its content with `taxonomyTree` and `ggtree`.

If you got stuck, you can look up chapters \@fref{datamanipulation}
and \@ref(fly-tree) on how to pick specific abundances and generate
row trees, respectively.


### Other elements

Try to extract some of the other TreeSE elements listed in chapter \@ref(containers).
However, such data are not always included.

1. Import the mia package, load peerj13075 with `data` and store it into a
   variable named `tse`.
2. Fetch the metadata of the TreeSE. Is there any iformation available?
3. Access the phylogenetic tree with `rowTree`. How big is it in terms of tips
   and nodes. If you like you can visualize it with `ggtree`.
4. Check if a sample tree is available with `colTree`, which is suitable for
   hierarchical or nested study designs.
5. If present, obtain the information on feature DNA sequences from the DNA
   sequence slot.


## Data manipulation


### Subsetting

1. Subset the TreeSE object to specific samples
2. Subset the TreeSE object to specific features
3. Subset the TreeSE object to specific samples and features



### Library sizes

1. Calculate library sizes
2. Subsample / rarify the counts (see: subsampleCounts)

Useful functions: nrow, ncol, dim, summary, table, quantile, unique, addPerCellQC, mergeFeaturesByRank

### Prevalent and core taxonomic features

1. Estimate prevalence for your chosen feature (row, taxonomic group)
2. Identify all prevalent features and subset the data accordingly 
3. Report the thresholds and the dimensionality of the data before and after subsetting
4. Visualize prevalence

Useful functions: getPrevalence, getPrevalentFeatures, subsetByPrevalentFeatures


### Data exploration

1. Summarize sample metadata variables. (How many age groups, how they are distributed? 0%, 25%, 50%, 75%, and 100% quantiles of library size?)
2. Create two histograms. Another shows the distribution of absolute counts, another shows how CLR transformed values are distributed.
3. Visualize how relative abundances are distributed between taxa in samples.

Useful functions: nrow, ncol, dim, summary, table, quantile, unique, transformAssay, ggplot, wilcox.test, mergeFeaturesByRank, plotAbundance

### Other functions

1. Merge data objects (merge, mergeSEs)
2. Melt the data for visualization purposes (meltAssay)


### Assay transformation

1. Import the mia package, load peerj13075 with `data` and store it into a
   variable named `tse`.
2. Transform the counts assay into relative abundances with `transformAssay` and
   store it into the TreeSE as an assay named `relabund` (see chapter \ref(assay-transform)).
3. Similarly, perform a clr transformation on the counts assay with a `pseudocount`
   of 1 and add it to the TreeSE as a new assay.
4. List the available assays by name with `assays`.
5. Access the clr assay and select a subset of its first 100 features and 10
   samples. Remember that assays are subsettable with `assay[row_idx, col_idx]`.
5. Take the same subset from the TreeSE, and check how this affects the individual
   transformed assays. TreeSE can also be subsetted with `tse[row_idx, col_idx]`.
6. **Extra**: If the data has phylogenetic tree, perform the phILR transformation.


## Abundance tables

### Taxonomic levels

1. Import the mia package, load peerj13075 with `data` and store it into a
   variable named `tse`.
2. List the available taxonomic ranks in the data with `taxonomyRanks`.
3. Agglomerate the data to Phylum level with `mergeFeaturesByRank` and the
   appropriate value for `Rank`.
4. Report the dimensions of the TreeSE before and after agglomerating. You can
   use `dim` for that.
5. **Extra**: Perform CLR transformation on the data. Does this affect agglomeration?
6. **Extra**: List full taxonomic information for a few selected taxa, such as
   `OTU1` and `OTU1368`. For that you can use `mapTaxonomy` on a specific subset
   of the TreeSE.


### Alternative experiments

1. Import the mia package, load GlobalPatterns with `data` and store it into a
   variable named `tse`.
2. Check the taxonomic ranks of the features with `taxonomyRanks`. What is the
   deepest taxonomic rank available?
3. Agglomerate the TreeSE to each taxonomic rank and store the resulting
   experiments as altExps. This can be performed automatically with `splitByRanks`.
4. Check the names of the generated altExps with `altExpNames` and retrieve a
   complete list with `altExps`.
5. Retrieve the data agglomerated by genus from the corresponding `altExp`.
   As for assays, you can access the desired altExp by name or index.
6. **Extra**: Split the data based on other features with `splitOn`.


## Community (alpha) diversity

### Estimation

1. Import the mia package, load GlobalPatterns with `data` and store it into a
   variable named `tse`.
2. Calculate multiple alpha diversity indices with `estimateDiversity` without
   any additional arguments.
3. Check the names of colData with `names`. Can you identify which columns
   contain the alpha diversity indices?
4. **Extra**: Agglomerate the TreeSE by phylum and compare the mean Shannon
   diversity of the original experiment with its agglomerated version. You can
   use `mergeFeaturesByRank` to perform agglomeration and `mean` to calculate the
   mean values of the respective columns in colData.

### Visualization 

1. Import the mia and scater packages, load GlobalPatterns with `data` and store
   it into a variable named `tse`.
2. Calculate Shannon diversity index and Faith's phylogenetic diversity with
   `estimateDiversity` and the appropriate arguments for `index`.
3. Make a boxplot of Shannon diversity on the y axis and sample type on the x
   axis with `plotColData`.
4. Repeat the previous point with Faith's phylogenetic diversity and compare the
   sample distributions of the two alpha diversity indices. How greatly do they
   differ?
5. **Extra**: Make a scatterplot of Shannon diversity on the y axis and Faith's
   phylogenetic diversity on the x axis with `plotColData`. Colour the points
   by sample type with the appropriate optional argument.
   
### Correlation

1. Import the mia and scater packages, load peerj13075 with `data` and store it
   into a variable named `tse`.
2. Calculate coverage and Shannon diversity index with `estimateDiversity` and
    the appropriate arguments for `index`.
4. Test the correlation between the two indices with `cor.test`. Remember that
   colData parameters are accessible with `tse$param_name`. Use Kendall tau
   coefficients as `method` to measure correlation. Is the correlation
   weak or strong, significant or not?
5. Make a scatterplot of Shannon diversity index on the y axis and coverage on
   the x axis. You can do that with `plotColData`. How do the two indices
   relate to one another?
6. **Extra**: Compute the library size of the samples by applying `colSums` to
   the counts assay of the TreeSE, and test the correlation of library size with
   Shannon diversity or coverage. Which index is more correlated with library size?
   
In this example, we inspected the correlation between two related variables,
also known as multicollinearity, and checked the correlation to library size,
which is part of quality control. However, the correlation between alpha diversity
and other numerical data about samples, such as participant's age and weight, also
represent an important analysis in several studies.

### Differences between groups

1. Import the mia package, load peerj13075 with `data` and store it into a
   variable named `tse`.
2. Calculate the Gini-Simpson diversity with `estimateDiversity` and the
   appropriate argument for `index`. Set `name` to `simpson`. You will use this
   name to access the diversity index from colData.
3. Inspect the Diet column in the colData. Determine how the samples are grouped
   in terms of diet. You can see the number of unique elements in a column with
   `unique`.
4. Test differences in Gini-Simpson diversity between different diets with
   `kruskal.test`. Remember that colData parameters are accessible with
   `tse$param_name`.
5. Is diversity significantly different between vegan and mixed diet? To
   visualize that, make a boxplot of Gini-Simpson diversity on the y axis and
   diet on the x axis with `plotColData`.
6. **Extra**: Repeat points 3 through 5, this time for age groups. Make sure
   that you are using an appropriate statistical test for the number of groups
   and distribution.


## Community similarity

### Reduced dimensions retrieval

1. Import the mia package, load enterotype with `data` and store it into a
   variable named `tse`.
2. List all available reduced dimensions with `reducedDims`. At this point, no
   reducedDims are likely found, because we haven't created any yet.
3. Perform PCA and store its output in the TreeSE by running
   `tse <- runPCA(tse, assay.type = "counts")`. Note that it is required to
   specify the assay on which dimensionality reduction should be conducted.
4. View the names of all available reduced dimensions with `reducedDimNames`.
   Has something new appeared?
5. **Extra**: Access the `PCA` reducedDim object with `reducedDim` and explore
   its content. How are the different dimensions stored? Try to extract an
   array with only the values from the second dimension by indexing the object
   with `[ , 2]`.

### Visualization basics with PCA

1. Import the mia and scater packages, load enterotype with `data` and store it
   into a variable named `tse`.
2. Perform a 3-component PCA based on the `counts` assay. You can use `runPCA`
   and set the optional arguments `ncomponents` and `assay.type` to the
   appropriate values.
3. Plot the first two dimensions of PCA with `plotReducedDim`, to which you
   should give the appropriate reducedDim name as the second argument. Note that
   by default only the first two dimensions are shown.
4. Check which information is stored in the ColData of the TreeSE. What would
   be worth visualizing in our coordination plot?
5. Make the same plot again, but this time colour the observations by Enterotype.
   You can do that by setting `colour_by` to the appropriate colname in the
   colData of the TreeSE.
6. **Extra**: Plot all three dimensions of PCA with `plotReducedDim` and the
   optional argument `ncomponents`. Colour observations by Enterotype. Which
   pair of dimensions best explains the variance between Enterotypes?


### Principal Coordinate Analysis (PCoA)

PCoA turns out to be particularly relevant for microbiome analysis, because
unlike PCA it can generate reduced dimensions from distances other than
Euclidean. There are several ecological distances to choose from and you can
find many of them under methods in the vignettes of `vegan::vegdist`.

1. Import the mia and scater packages, load enterotype with `data` and store it
   into a variable named `tse`.
2. Transform the counts assay to relative abundances with `transformAssay`.
3. Perform a Multi-Dimensional Scaling (MDS) based on the relative abundance
   assay in terms of Bray-Curtis dissimilarity. You can use `runMDS` with the
   compulsory argument `FUN = vegan::vegdist`.
4. Plot the first two dimensions of PCA with `plotReducedDim`, to which you
   should give the appropriate reducedDim name as the second argument. Colour
   the observations by Enterotype with `colour_by`.
5. **Extra**: Perform MDS again with `runMDS`, but this time use Jaccard
   dissimilarity. The distance metric to use can be defined with the optional
   argument `method`, choosing from the methods in `?vegan::vegdist`. If you
   don't want to overwrite the reducedDim object made in point 3, set `name` to
   a name of your choice. Visualize and compare it to the plot from point 4.
   
Good job! You are now able to produce and visualize reduced dimensions of a
TreeSE. `runMDS` is actually one of several algorithms for PCoA and dimensionality
reduction, which you can find in section \@ref(other-ord-methods).

### PERMANOVA analysis

In this exercise we focus on studying the weight of variables on the microbiome
composition. Significance of each variable on beta diversity is tested with
PERMANOVA (point 4) and the homogeneity assumption is also be controlled with a
PERMDISP2 analysis (point 5).

1. Import the mia and vegan packages, load peerj13075 with `data` and store it
   into a variable named `tse`.
2. Transform the `counts` assay into relative abundances with `transformAssay`.
3. Extract the relative abundance assay, transpose it and save it into
   a variable named `relabund_assay`.
4. Perform PERMANOVA with `adonis2` to see how much Diet can explain the relative
   abundance assay (`formula = relabund_assay ~ Diet`) in terms of Bray-Curtis
   dissimilarity (`method = "bray"`). Also set `data = colData(tse)`
   `by = "margin"` and `permutations = 99`. What do the results tell you about
   Diet with respect to beta diversity?
5. **Extra**: Test homogeneity of distribution across Diet groups with
   `anova(betadisper(my_mat), my_groups`, where `my_mat` is the Bray-Curtis
   dissimilarity matrix of `relabund_assay` calculated with
   `vegdist(relabund_assay, "bray")` and `my_groups` is the vector of Diet values
   obtained from the colData of the TreeSE.
   
Well done! You went through testing the effect and significance of Diet on beta
diversity. Keep in mind that the formula fed to `adonis2` can take more than one
independent variable, so that you can also (and very often should) include
covariates of your studies.

### Redundancy analysis (RDA)

Here we apply RDA, an ordination method that provides dimensions with
the largest variation between the data based in the light of the specified
variables (point 3). The results of RDA are usually assessed with PERMANOVA
(point 5) and the homogeneity assumption should be checked as in the previous
exercise. This is a relatively complex procedure, but the way this is broken
down into steps below will hopefully make more sense.

1. Import the mia and vegan packages, load peerj13075 with `data` and store it
   into a variable named `tse`.
2. Transform the `counts` assay into relative abundances with `transformAssay`.
3. Perform RDA with `calculateRDA` to see how much Diet can explain the relative
   abundance assay (`formula = assay ~ Diet` and `assay.type = relabundance`) in
   terms of Bray-Curtis dissimilarity (`method = "bray"`).
4. Extract the RDA dimensions from the appropriate reducedDim slot with
   `attr(reducedDim(tse, "RDA"), "rda)` and store it into `rda`.
5. Test the effect and significance of Diet on beta diversity by PERMANOVA with
   `anova.cca`. Feed this function with `rda` and set `by = "margin"` and
   `permutations = 99`, respectively. What do the results tell you about Diet?
6. **Extra**: Check what other parameters are stored in the colData of peerj13075,
   add them to the formula (`formula = assay ~ Diet + ...`) of `calculateRDA`
   and proceed to see how that changes the results of PERMANOVA.

Well done! You went through an RDA analysis followed by significance testing
with PERMANOVA and BETADISPER2. In the next exercise we'll go deeper quantify
the contributions to beta diversity.

### Beta diversity analysis

This exercise prompts you to implement a workflow with distance-based RDA
(dbRDA). You can refer to chapter \@ref(dbrda-workflow) for a step-by-step
walkthrough, which may be simplified in the future.

1. Import the mia and vegan packages, load peerj13075 with `data` and store it
   into a variable named `tse`.
4. Create dbRDA with Bray-Curtis dissimilarities on relative abundances. Use PERMANOVA. Can differences between samples be explained with variables of sample meta data? 
5. Analyze diets' association on beta diversity. Calculate dbRDA and then PERMANOVA. Visualize coefficients. Which taxa's abundances differ the most between samples? 
6. Interpret your results. Is there association between community composition and location? What are those taxa that differ the most; find information from literature.

Useful functions: runMDS, runRDA, anova.cca, transformAssay, mergeFeaturesByRank, ggplot, plotReducedDim, vegan::adonis2




## Differential abundance

### Univariate analyses

1. Get the abundances for an individual feature (taxonomic group / row)
2. Visualize the abundances per group with boxplot / jitterplot
3. Is the difference significant (Wilcoxon test)?
4. Is the difference significant (linear model with covariates)? 
5. How do transformations affect the outcome (log10, clr..)?
6. Get p-values for all features (taxa), for instance with a for loop
7. Do multiple testing correction
9. Compare the results from different tests with a scatterplot

Useful functions: [], ggplot2::geom_boxplot, ggplot2::geom_jitter, wilcox.test, lm.test, transformAssay, p.adjust


### Differential abundance analysis

1. install the latest development version of mia from GitHub.
2. Load experimental dataset from mia.
3. Compare abundances of each taxa between groups. First, use Wilcoxon or Kruskall-Wallis test. Then use some other method dedicated to microbiome data. 
4. Summarize findings by plotting a taxa vs samples heatmap. Add column annotation that tells the group of each sample, and row annotation that tells whether the difference of certain taxa was statistically significant.
5. Choose statistically significant taxa and visualize their abundances with boxplot & jitterplot.


Useful functions: wilcox.test, kruskal.test, ggplot, pheatmap, ComplexHeatMap::Heatmap, ancombc, aldex2, maaslin2, mergeFeaturesByRank, transformAssay, subsetByPrevalentFeatures



## Visualization

### Multivariate ordination

1. Load experimental dataset from mia.
2. Create PCoA with Bray-Curtis dissimilarities
3. Create PCA with Aitchison dissimilarities
4. Visualize and compare both
5. Test other transformations, dissimilarities, and ordination methods

Useful functions: runMDS, runNMDS, transformAssay, ggplot, plotReducedDim


### Heatmap visualization

1. Load experimental dataset from mia.
2. Visualize abundances with heatmap
3. Visualize abundances with heatmap after CLR + Z transformation 

See the OMA book for examples.



## Multiomics

### Basic exploration

Here we learn how to conduct preliminary exploration on a MAE, using
HintikkaXOData as an example dataset.

1. Import the mia package, load HintikkaXOData with `data` and store it into a
   variable named `mae`.
2. Which experiments make up the MAE? How many samples and features are contained
   by each experiment? You can get a summary for all experiments with `experiments`,
   and check for each individual experiment with `dim`, `nrow` and `ncol`.
3. What are the names of the features and samples of the different experiments?
   You can see that with `rownames` and `colnames`, respectively.
4. What information is known about the samples? Remember that information about
   samples is stored in the `colData` of the MAE.
5. **Extra**: How do the samples of the individual experiments map to the
   columns of the MAE? You can find the sample mapping in the `sampleMap`
   of the MAE.

So far so good. You explored a MAE and its experiments, getting a taste of how
information is organized in its internal structure.

### Experiment agglomeration

Here we learn how to manipulate an experiment contained by a MAE and save the
new modified version of the experiment in a suitable place (the altExp slot).

1. Import the mia package, load HintikkaXOData with `data` and store it into a
   variable named `mae`.
2. Agglomerate the microbiota experiment by Genus and store the output into the
   `altExp` slot of the microbiota experiment, with the custom name
   `microbiota_genus`.
3. How many features remain after agglomerating? What are their names?
4. **Extra**: create one more alternative experiment named
   `prevalent_microbiota_family`, which contains the microbiota experiment
   agglomerated by Family with a prevalence threshold of 10%. You can
   agglomerate and in parallel select by prevalence with
   `mergeFeaturesByPrevalence`.

Good job! You agglomerated one of the experiments in the MAE and stored it as
an alternative experiment.

### Experiment transformation

We proceed with an exercise on a different type of data manipulation, that is,
transformation of assays of individual experiments in the MAE.

1. Import the mia package, load HintikkaXOData with `data` and store it into a
   variable named `mae`.
2. What assays are contained by each individual experiment? You can check their
   names with `assays`.
3. Apply a log10 transformation to the assay of the metabolite experiment.
   For that you can use `transformAssay` and don't forget to specify the assay
   to be transformed with the argument `assay.type`.
4. Apply a CLR transformation to the counts assay of the microbiota experiment.
   To ensure non-null values in the assay, set `pseudocount` equal to 1.
   
You made it! You learnt how to apply different transformations to the assays
of individual experiments in a MAE with `transformAssay`, specifying optional
arguments based on the used method.

### Assay extraction

The following exercise walks you through disassembling a MAE object in
order to retrieve a specific assay, or to store its components as
multiple separate csv files.

1. Import the mia package, load HintikkaXOData with `data` and store it into a
   variable named `mae`.
2. Extract the individual metabolite experiment from the MAE into a distinct
   TreeSE object named `metabolites`.
3. Which and how many assays are contained by `metabolites`? You can check that
   with `assays` or `assayNames`.
4. Write a csv file for the nmr assay with `write.csv`. You can access an
   individual assay of a TreeSE with `assay` by specifying the name of the
   desired assay.
5. **Extra**: Repeat step 1 thorugh 4 also for the microbiota and biomarkers
   experiments, so that a completely disassembled version of the MAE is available.
6. **Extra**: Besides experiments, MAEs also include a sampleData and a sampleMap,
   which are accessible with `colData(mae)` and `sampleMap(mae)`, respectively.
   Save also each of these two elements into a csv file.

Well done! You just splitted a MAE into its components and stored them as csv files.
[This script](https://github.com/JuliaTurkuDataScience/MicrobiomeAnalysis.jl/blob/main/src/assets/XO_preprocess.R)
shows a possible approach.

### MAE reconstruction

Next, we will try to reconstruct the same MAE from the files you created.
Make sure you know their names and location! Alternatively, you can fetch or
download the CSV files in
[this directory](https://github.com/microbiome/data/tree/main/HintikkaXOData)
with the readily disassembled components of HintikkaXOData.

1. Read in the csv files containing assays with `read.csv` and save each
   of them into a variable named `<assay name>_assays`.
2. Create one TreeSE from each assays object with the `TreeSummarizedExperiment`
   function, as explained in [this exercise](#construct-TreeSE).
3. Read in the sampleData and the sampleMap and store them into the
   variables `sample_data` and `sample_map`, respectively.
4. Combine the components with `MultiAssayExperiment`, where the first argument
   is an `ExperimentList` (for now include only the microbiota and metabolites
   TreeSEs), the second is colData and the third is sampleMap.
5. Make sure that the MAE experiments are identical to the original TreeSEs. You
   can do that qualitatively by checking their `head` and quantitatively by 
   looking at their `dim`.
6. **Extra**: Add the biomarkers TreeSE as a new experiment to the MAE.
   Note that new experiments can be added to a MAE through simple concatenation
   with `c(mae, experiment)`.

Good job! Now you are aware of how MAEs are built and we can proceed to some
analytical exercises.

### Cross-correlation analysis

Now we will perform a cross-correlation analysis between two of the
experiments in the MAE.

1. Import the mia package, load HintikkaXOData with `data` and store it into a
   variable named `mae`.
2. Analyze correlations between the microbiota and the biomarkers experiments
   with `getExperimentCrossAssociation`. Don't forget to specify the experiments
   you want to compare with the arguments `experiment1` and `experiment2`, and
   which of their assays with `assay.type1` and `assay.type2`.
3. What does the output look like? By default, correlation is measured in terms
   of Kendall tau coefficients. Repeat point 2, but this time change `method`
   to Spearman coefficients.
4. Are you able to infer significance from the output? In order to also obtain
   p-values from the cross-correlation analysis, repeat point 2 with the
   additional argument `test_significance = TRUE`.
5. Visualize results with a heatmap similarly to the example in section
   \@ref(cross-correlation). Do you see any significant correlations?
   Interpret your results.
6. **Extra**: Perform cross-correlation analysis between the remaining
   experiments (microbiota vs metabolites and metabolites vs biomarkers) and
   visualize results with heatmaps.
   
Great job! You performed a cross-correlation analysis between two experiments of
a MAE and visualized the results with a heatmap. You are also able to customise
the correlation method and significance testing used for the analysis.
