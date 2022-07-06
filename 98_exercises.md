# Exercises {#exercises}

Here you can find assignments on different topics. 

TIPS:

   - Add comments that explain what each line or lines of code do. This helps you and others to understand your code and find bugs. Furthermore, it is easier for you to reuse the code, and it promotes transparency.
   - Interpret results by comparing them to literature. List main findings, so that results can easily be understood by others without advanced data analytics knowledge.
   - Avoid hard-coding. Use variables which get values in the beginning of the pipeline. That way it is easier for you to modify parameters and reuse the code.

## Workflows

### Reproducible reporting

1. Create an R Markdown file.
2. Add a code chunk and give a name for it.
3. Import e.g., iris dataset, and create a dotplot with a title.
4. Create another code chunk and plot. 
5. Adjust the size of the figure, and choose from the options that the code chunk will not be shown in the report.
6. Add some text.
7. Add an R command within the text.
8. Create an HTML file from the Rmd file.

## Data containers

### Introduction to TreeSummarizedExperiment (TreeSE)

Import data from CSV files to TreeSE.

1. Import data.
2. Construct a TreeSE.
3. Check that importing is done correctly. E.g., you can choose random samples and features,
and check that their values equal between raw files and TreeSE.

Useful functions: DataFrame, TreeSummarizedExperiment, matrix, rownames, colnames, SimpleList

### TreeSummarizedExperiment Data exploration

1. install the latest development version of mia from GitHub.
2. Load experimental dataset from mia.
3. What are the dimensions? (How many samples there are, and how many taxa in each taxonomic rank)?
4. Calculate library size.
5. Summarize sample metadata variables. (How many age groups, how they are distributed? 0%, 25%, 50%, 75%, and 100% quantiles of library size?)
6. Create two histograms. Another shows the distribution of absolute counts, another shows how CLR transformed values are distributed.
7. Visualize how relative abundances are distributed between taxa in samples.


Useful functions: nrow, ncol, dim, summary, table, quantile, unique, transformSamples, ggplot, wilcox.test, addPerCellQC, agglomerateByRank, plotAbundance



### Introduction to MultiAssayExperient (MAE)

1. Create TreeSE data containers from individual CSV files.
2. Combine TreeSE into MAE.
3. Check that each individual experiment of MAE equals corresponding TreeSE.
4. Take a subset of MAE (e.g., 10 first samples), and observe the subsetted MAE.

Useful functions: DataFrame, TreeSummarizedExperiment, matrix, rownames, colnames, MultiAssayExperiment, ExperimentList, SimpleList



## Data manipulation

### Prevalent and core features

1. Estimate prevalence for your chosen feature (row, taxonomic group)
2. Identify all prevalent features and subset the data accordingly 
3. Report the thresholds and the dimensionality of the data before and after subsetting
4. Visualize prevalence

Useful functions: getPrevalence, getPrevalentTaxa


### Taxonomic levels

1. List the available taxonomic ranks in the data
2. Split the data into alternative experiments by rank (altExp)
3. Merge the data to Phylum level; report dimensionality before and after this

Useful functions: taxonomyRanks, agglomerateByRank, mergeRows, splitByRank


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


## Alpha diversity

### Alpha diversity basics

1. Calculate alpha diversity indices
2. Test if data agglomeration to higher taxonomic ranks affects the indices
3. Look for differences in alpha diversity between groups or correlation with a continuous variable


### Alpha diversity 

1. Estimate Shannon diversity for the data
2. Try also another diversity index and compare the results with a scatterplot
3. Compare Shannon diversity between groups (boxplot)
4. Is diversity significantly different between vegan and mixed diet?
5. Calculate and visualize library size, compare with diversity

Useful functions: estimateDiversity, colSums, geom_point, geom_boxplot





## Community composition

### Beta diversity

1. Visualize community variation with different methods (PCA, MDS, NMDS, etc.) with plotReduceDim and with different dissimilarities and transformations,plot also other than the first two axes.
2. Use PERMANOVA to test differences in beta diversity. You can also try including continuous and/or categorical covariates
3. If there are statistical differences in PERMANOVA, test PERMDISP2 (betadisper function)
4. Do DMM clustering
5. Try RDA to test the variance explained by external variables



### Basics of community composition

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


## Multiomics

### Introduction to multiomics

1. Load experimental dataset from microbiomeDataSets (e.g., HintikkaXOData).
2. Analyze correlations between experiments. (Taxa vs lipids, Taxa vs biomarkers, Lipids vs biomarkers)
3. Agglomerate taxa data.
4. Apply CLR to taxa data, apply log10 to lipids and biomarkers.
5. Perform cross-correlation analyses and visualize results with heatmaps. (Use Spearman coefficients)
6. Is there significant correlations? Interpret your results.

Useful functions: pheatmap, ComplexHeatMap::Heatmap, ggplot, transformSamples, testExperimentCrossAssociation






