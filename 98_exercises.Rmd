# Exercises {#exercises}

Here you can find assignments on different topics. 

**Tips for exercises:**

   - Add comments that explain what each line or lines of code do. This helps you and others to understand your code and find bugs. Furthermore, it is easier for you to reuse the code, and it promotes transparency.
   - Interpret results by comparing them to literature. List main findings, so that results can easily be understood by others without advanced data analytics knowledge.
   - Avoid hard-coding. Use variables which get values in the beginning of the pipeline. That way it is easier for you to modify parameters and reuse the code.

## Workflows

### Reproducible reporting

1. Create a new Quarto file

* Rstudio has ready-made templates for this
* [Creating Quarto Documents](https://quarto.org/docs/tools/rstudio.html#creating-documents)

2. Add a code chunk and name it.
3. Render (or knit) the file into pdf or html format (Hint: [Quarto Rstudio rendering](https://quarto.org/docs/tools/rstudio.html#render-and-preview))
4. Import e.g., iris dataset, and add a dotplot with a title (Hint: [geom_dotplot](https://ggplot2.tidyverse.org/reference/geom_dotplot.html))
5. Create another code chunk and plot. 
6. Adjust figure size and hide the code chunk from output report.

  * [Sizing](https://quarto.org/docs/authoring/diagrams.html#sizing)
  * [chunk-options](https://quarto.org/docs/computations/r.html#chunk-options)

7. Add some text.
8. Add R commands within the text (Hint: [code-chunks](https://quarto.org/docs/visual-editor/technical.html#code-chunks))
9. Update HTML file from the qmd file (Hint: [Quarto Rstudio rendering](https://quarto.org/docs/tools/rstudio.html#render-and-preview))


For tips on Quarto, see [Quarto tutorial](https://quarto.org/docs/authoring/markdown-basics.html)


## TreeSummarizedExperiment: basic components

You can also check the [function reference in the mia package](https://microbiome.github.io/mia/reference/index.html)

### Basic summaries

1. Load experimental dataset from mia (e.g. `peerj13075` with the `data()` command; see OMA [section 2.4 Demonstration Data](https://microbiome.github.io/OMA/containers.html#example-data); also see [loading experimental data](https://microbiome.github.io/OMA/containers.html#assay-data)).
2. Check a summary about the TreeSE object loaded (Hint: `summary()`)
3. What are the dimensions? (How many samples there are, and how many taxa in each taxonomic rank?) (Hint: material in [Section 2](https://microbiome.github.io/OMA/containers.html#data-science-framework) may help)
4. List sample and features names for the data (rownames, colnames..)

### Taxonomic abundance table (assay)

1. (Load example data)
2. Fetch the list of available assays (Hints: [assayNames](https://microbiome.github.io/OMA/containers.html#assay-data))
3. Fetch the `counts` assay, and show part of it. (Hint: [assay-data](https://microbiome.github.io/OMA/containers.html#assay-data))


### Sample side information 

1. (Load example data)
2. Fetch and show data about samples (Hint: [colData](https://microbiome.github.io/OMA/containers.html#coldata))
3. Get abundance data for all taxa for a specific sample (sample names: function `colnames(tse)`)

  * [example](https://microbiome.github.io/OMA/taxonomic-information.html#abundances-of-all-taxa-in-specific-sample)


### Feature side information

1. (Load example data)
2. Fetch and show data on the (taxonomic) features of the analyzed samples (Hint: [rowData](https://microbiome.github.io/OMA/containers.html#rowdata))
3. Get abundance data for all samples given a specific features (Hint: [example](https://microbiome.github.io/OMA/taxonomic-information.html#abundances-of-specific-taxa-in-all-samples))

Optional:

4. Create taxonomy tree based on the taxonomy mappings display its information:

 * [generate a taxonomic tree on the fly](https://microbiome.github.io/OMA/taxonomic-information.html#generate-a-taxonomic-tree-on-the-fly)
 * [rowtree](https://microbiome.github.io/OMA/containers.html#rowtree)


### Other components

Try to extract some of the [other TreeSE elements](https://f1000research.com/articles/9-1246/v2). These are not always included:

* Experiment metadata
* Sample tree (colTree)
* Phylogenetic tree (rowTree)
* Feature sequences information (DNA sequence slot)


## TreeSummarizedExperiment: data import

Import data from CSV files to TreeSE (see shared data folder for example data sets).

1. Import the data files in R
2. Construct a TreeSE data object (see [Ch. 2](https://microbiome.github.io/OMA/containers.html#loading-experimental-microbiome-data))
3. Check that importing is done correctly. E.g., choose random samples and features,
and check that their values equal between raw files and TreeSE.

Useful functions: DataFrame, TreeSummarizedExperiment, matrix, rownames, colnames, SimpleList

Optional:

4. Import data from another format (functions: loadFromMetaphlan | loadFromMothur | loadFromQIIME2 | makeTreeSummarizedExperimentFromBiom | makeTreeSummarizedExperimentFromDADA2 ...)
5. Try out conversions between TreeSE and phyloseq data containers (makeTreeSummarizedExperimentFromPhyloseq; makephyloseqFromTreeSummarizedExperiment)





## Data manipulation


### Subsetting

1. Subset the TreeSE object to specific samples
2. Subset the TreeSE object to specific features
3. Subset the TreeSE object to specific samples and features



### Library sizes

1. Calculate library sizes
2. Subsample / rarify the counts (see: subsampleCounts)

Useful functions: nrow, ncol, dim, summary, table, quantile, unique, addPerCellQC, agglomerateByRank

### Prevalent and core taxonomic features

1. Estimate prevalence for your chosen feature (row, taxonomic group)
2. Identify all prevalent features and subset the data accordingly 
3. Report the thresholds and the dimensionality of the data before and after subsetting
4. Visualize prevalence

Useful functions: getPrevalence, getPrevalentTaxa, subsetByPrevalentTaxa


### Explore the data

1. Summarize sample metadata variables. (How many age groups, how they are distributed? 0%, 25%, 50%, 75%, and 100% quantiles of library size?)
2. Create two histograms. Another shows the distribution of absolute counts, another shows how CLR transformed values are distributed.
3. Visualize how relative abundances are distributed between taxa in samples.

Useful functions: nrow, ncol, dim, summary, table, quantile, unique, transformSamples, ggplot, wilcox.test, agglomerateByRank, plotAbundance

### Other functions

1. Merge data objects (merge, mergeSEs)
2. Melt the data for visualization purposes (meltAssay)


## Transformations

### Transformations

1. Transform abundance data with relative abundance and add a relative abundance assay (see [data-transformation](https://microbiome.github.io/OMA/taxonomic-information.html#data-transformation))
2. Transform abundance data with clr transformation and add a new assay
3. List the available assays by name
4. Pick one of the assays and show a subset of it
5. Subset the entire TreeSE data object, and check how this affects individual (transformed) assays

Optional:

6. If the data has phylogenetic tree, perform the phILR transformation


## Taxonomic levels

### Taxonomic levels

1. List the available taxonomic ranks in the data
2. Merge the data to Phylum level
3. Report dimensionality before and after aggomeration

Optional:

4. Perform CLR transformation on the data; does this affect agglomeration?
5. List full taxonomic information for some given taxa (Hint: [mapTaxonomy](https://microbiome.github.io/mia/reference/taxonomy-methods.html))

Useful functions: [taxonomyRanks](https://microbiome.github.io/mia/reference/taxonomy-methods.html), agglomerateByRank, mergeRows


### Alternative experiments (altExp)

1. Create taxonomic abundance tables for all different levels (splitByRanks)
2. Check the available alternative experiment (altExp) names before and after splitByRanks
3. Pick specific "experiment" (taxonomic rank) from specific altExp; and a specific assay

Optional:

4. Split the data based on other features (splitOn)



## Alpha diversity

### Alpha diversity basics

1. Calculate alpha diversity indices
2. Test if data agglomeration to higher taxonomic ranks affects the indices
3. Look for differences in alpha diversity between groups or correlation with a continuous variable

Useful functions: estimateDiversity, colSums, agglomerateByRank, kruskal.test, cor


### Alpha diversity extra

1. Estimate Shannon diversity for the data
2. Try also another diversity index and compare the results with a scatterplot
3. Compare Shannon diversity between groups (boxplot)
4. Is diversity significantly different between vegan and mixed diet?
5. Calculate and visualize library size, compare with diversity

Useful functions: estimateDiversity, colSums, geom_point, geom_boxplot




## Community composition

### Beta diversity basics

1. Visualize community variation with different methods (PCA, MDS, NMDS, etc.) with plotReduceDim and with different dissimilarities and transformations,plot also other than the first two axes.
2. Use PERMANOVA to test differences in beta diversity. You can also try including continuous and/or categorical covariates
3. If there are statistical differences in PERMANOVA, test PERMDISP2 (betadisper function)
4. Do clustering
5. Try RDA to test the variance explained by external variables


### Beta diversity extra

1. Install the latest development version of mia from GitHub.
2. Load experimental dataset from mia.
3. Create a PCoA with Aitchison dissimilarities. How much coordinate 1 explains the differences? How about coordinate 2?
4. Create dbRDA with Bray-Curtis dissimilarities on relative abundances. Use PERMANOVA. Can differences between samples be explained with variables of sample meta data? 
5. Analyze diets' association on beta diversity. Calculate dbRDA and then PERMANOVA. Visualize coefficients. Which taxa's abundances differ the most between samples? 
6. Interpret your results. Is there association between community composition and location? What are those taxa that differ the most; find information from literature.

Useful functions: runMDS, runRDA, anova.cca, transformSamples, agglomerateByRank, ggplot, plotReducedDim, vegan::adonis2







## Visualization

### Multivariate ordination

1. Load experimental dataset from mia.
2. Create PCoA with Bray-Curtis dissimilarities
3. Create PCA with Aitchison dissimilarities
4. Visualize and compare both
5. Test other transformations, dissimilarities, and ordination methods

Useful functions: runMDS, runNMDS, transformSamples, ggplot, plotReducedDim


### Heatmap visualization

1. Load experimental dataset from mia.
2. Visualize abundances with heatmap
3. Visualize abundances with heatmap after CLR + Z transformation 

See the OMA book for examples.


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

Useful functions: [], ggplot2::geom_boxplot, ggplot2::geom_jitter, wilcox.test, lm.test, transformSamples, p.adjust


### Differential abundance analysis

1. install the latest development version of mia from GitHub.
2. Load experimental dataset from mia.
3. Compare abundances of each taxa between groups. First, use Wilcoxon or Kruskall-Wallis test. Then use some other method dedicated to microbiome data. 
4. Summarize findings by plotting a taxa vs samples heatmap. Add column annotation that tells the group of each sample, and row annotation that tells whether the difference of certain taxa was statistically significant.
5. Choose statistically significant taxa and visualize their abundances with boxplot & jitterplot.

Useful functions: wilcox.test, kruskal.test, ggplot, pheatmap, ComplexHeatMap::Heatmap, ancombc, aldex2, maaslin2, agglomerateByRank, transformSamples, transformFeatures, subsetByPrevalentTaxa




## Multiomics


### Introduction to MultiAssayExperiment (MAE)

1. Create TreeSE data containers from individual CSV files.
2. Combine TreeSE into MAE.
3. Check that each individual experiment of MAE equals corresponding TreeSE.
4. Take a subset of MAE (e.g., 10 first samples), and observe the subsetted MAE.

Useful functions: DataFrame, TreeSummarizedExperiment, matrix, rownames, colnames, MultiAssayExperiment, ExperimentList, SimpleList


### Introduction to multiomics

1. Load experimental dataset from microbiomeDataSets (e.g., HintikkaXOData).
2. Analyze correlations between experiments. (Taxa vs lipids, Taxa vs biomarkers, Lipids vs biomarkers)
3. Agglomerate taxa data.
4. Apply CLR to taxa data, apply log10 to lipids and biomarkers.
5. Perform cross-correlation analyses and visualize results with heatmaps. (Use Spearman coefficients)
6. Is there significant correlations? Interpret your results.

Useful functions: pheatmap, ComplexHeatMap::Heatmap, ggplot, transformSamples, testExperimentCrossAssociation






