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




# Community similarity {#community-similarity}

Where alpha diversity focuses on community variation within a
community (sample), beta diversity quantifies (dis-)similarites
between communities (samples). Some of the most popular beta diversity
measures in microbiome research include Bray-Curtis index (for
compositional data), Jaccard index (for presence / absence data,
ignoring abundance information), Aitchison distance (Euclidean
distance for clr transformed abundances, aiming to avoid the
compositionality bias), and the Unifrac distances (that take into
account the phylogenetic tree information). Only some of the commonly
used beta diversity measures are actual _distances_; this is a
mathematically well-defined concept and many ecological beta diversity
measures, such as Bray-Curtis index, are not proper distances.
Therefore, the term dissimilarity or beta diversity is commonly used.

Technically, beta diversities are usually represented as `dist`
objects, which contain triangular data describing the distance between
each pair of samples. These distances can be further subjected to
ordination. Ordination is a common concept in ecology that aims to
reduce the dimensionality of the data for further evaluation or
visualization. Ordination techniques aim to capture as much of
essential information in the data as possible in a lower dimensional
representation.  Dimension reduction is bound to loose information but
the common ordination techniques aim to preserve relevant information
of sample similarities in an optimal way, which is defined in
different ways by different methods. [TODO add references and/or link
to ordination chapter instead?]

Some of the most common ordination methods in microbiome research
include Principal Component Analysis (PCA), metric and non-metric
multi-dimensional scaling (MDS, NMDS), The MDS methods are also known
as Principal Coordinates Analysis (PCoA). Other recently popular
techniques include t-SNE and UMAP. 


## Explained variance

The percentage of explained variance is typically shown for PCA
ordination plots. This quantifies the proportion of overall variance
in the data that is captured by the PCA axes, or how well the
ordination axes reflect the original distances.

Sometimes a similar measure is shown for MDS/PCoA. The interpretation
is generally different, however, and hence we do not recommend using
it. PCA is a special case of PCoA with Euclidean distances.  With
non-Euclidean dissimilarities PCoA uses a trick where the pointwise
dissimilarities are first cast into similarities in a Euclidean space
(with some information loss i.e. stress) and then projected to the
maximal variance axes. In this case, the maximal variance axes do not
directly reflect the correspondence of the projected distances and
original distances, as they do for PCA.

In typical use cases, we would like to know how well the ordination
reflects the original similarity structures; then the quantity of
interest is the so-called "stress" function, which measures the
difference in pairwise similarities between the data points in the
original (high-dimensional) vs. projected (low-dimensional) space.

Hence, we propose that for PCoA and other ordination methods, users
would report relative stress (varies in the unit interval; the smaller
the better). This can be calculated as shown below. For further
examples, check the [note from Huber
lab](https://www.huber.embl.de/users/klaus/Teaching/statisticalMethods-lab.pdf).



```r
# Example data
library(mia)
data(GlobalPatterns, package="mia")

# Data matrix (features x samples)
tse <- GlobalPatterns
tse <- transformCounts(tse, method = "relabundance")

# Add group information Feces yes/no
colData(tse)$Group <- colData(tse)$SampleType=="Feces"

# Quantify dissimilarities in the original feature space
library(vegan)
x <- assay(tse, "relabundance") # Pick relabunance assay separately
d0 <- as.matrix(vegdist(t(x), "bray"))

# PCoA Ordination
pcoa <- as.data.frame(cmdscale(d0, k = 2))
names(pcoa) <- c("PCoA1", "PCoA2")

# Quantify dissimilarities in the ordination space
dp <- as.matrix(dist(pcoa))

# Calculate stress i.e. relative difference in the original and
# projected dissimilarities
stress <- sum((dp - d0)^2)/sum(d0^2)
```


Shepard plot visualizes the original versus projected (ordination)
dissimilarities between the data points:


```r
ord <- order(as.vector(d0))
df <- data.frame(d0 = as.vector(d0)[ord],
                  dmds = as.vector(dp)[ord])

library(ggplot2)
ggplot(aes(x = d0, y = dmds), data=df) + 
       geom_smooth() +
       geom_point() +       
       labs(title = "Shepard plot",
       x = "Original distance",
       y = "MDS distance",       
            subtitle = paste("Stress:", round(stress, 2))) +
  theme_bw()
```

![](20_beta_diversity_files/figure-latex/shepard-1.png)<!-- --> 


## Community comparisons by beta diversity analysis

A typical comparison of community composition starts with a visual
comparison of the groups on a 2D ordination.

Then we estimate relative abundances and MDS ordination based on
Bray-Curtis (BC) dissimilarity between the groups, and visualize the
results.

In the following examples dissimilarities are calculated by 
functions supplied to the `FUN` argument. This function can be defined by
the user. It must return a `dist` function, which can then be used to
calculate reduced dimensions either via ordination methods (such as MDS
or NMDS), and the results can be stored in the `reducedDim`.

This entire process is wrapped in the `runMDS` and `runNMDS`
functions.


```r
library(scater)

# Bray-Curtis is usually applied to relative abundances
tse <- transformCounts(tse, method = "relabundance")
# Perform PCoA
tse <- runMDS(tse, FUN = vegan::vegdist, method = "bray", name = "PCoA_BC", assay.type = "relabundance")
```

Sample similarities can be visualized on a lower-dimensional display
(typically 2D) using the `plotReducedDim` function in the `scater`
package. This provides also further tools to incorporate additional
information using variations in color, shape or size. Are there
visible differences between the groups?


```r
# Create ggplot object
p <- plotReducedDim(tse, "PCoA_BC", colour_by = "Group")

# Add explained variance for each axis
e <- attr(reducedDim(tse, "PCoA_BC"), "eig");
rel_eig <- e/sum(e[e>0])		  
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""))

print(p)
```

![(\#fig:plot-mds-bray-curtis)MDS plot based on the Bray-Curtis distances on the GlobalPattern dataset.](20_beta_diversity_files/figure-latex/plot-mds-bray-curtis-1.png) 




With additional tools from the `ggplot2` universe, comparisons can be 
performed informing on the applicability to visualize sample similarities in a 
meaningful way.


```r
tse <- runMDS(tse, FUN = vegan::vegdist, name = "MDS_euclidean",
             method = "euclidean", assay.type = "counts")
tse <- runNMDS(tse, FUN = vegan::vegdist, name = "NMDS_BC")
```

```
## initial  value 47.733208 
## iter   5 value 33.853364
## iter  10 value 32.891200
## final  value 32.823570 
## converged
```

```r
tse <- runNMDS(tse, FUN = vegan::vegdist, name = "NMDS_euclidean",
               method = "euclidean")
```

```
## initial  value 31.882673 
## final  value 31.882673 
## converged
```

```r
plots <- lapply(c("PCoA_BC", "MDS_euclidean", "NMDS_BC", "NMDS_euclidean"),
                plotReducedDim,
                object = tse,
                colour_by = "Group")

library(patchwork)
plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] +
  plot_layout(guides = "collect")
```

![(\#fig:plot-mds-nmds-comparison)Comparison of MDS and NMDS plots based on the Bray-Curtis or euclidean distances on the GlobalPattern dataset.](20_beta_diversity_files/figure-latex/plot-mds-nmds-comparison-1.png) 

The _Unifrac_ method is a special case, as it requires data on the
relationship of features in form on a `phylo` tree. `calculateUnifrac`
performs the calculation to return a `dist` object, which can again be
used within `runMDS`.



```r
library(scater)
tse <- runMDS(tse, FUN = mia::calculateUnifrac, name = "Unifrac",
              tree = rowTree(tse),
              ntop = nrow(tse),
             assay.type = "counts")
```


```r
plotReducedDim(tse, "Unifrac", colour_by = "Group")
```

![(\#fig:plot-unifrac)Unifrac distances scaled by MDS of the GlobalPattern dataset.](20_beta_diversity_files/figure-latex/plot-unifrac-1.png) 

## Other ordination methods

Other dimension reduction methods, such as `PCA`, `t-SNE` and `UMAP` are 
inherited directly from the `scater` package.


```r
tse <- runPCA(tse, name = "PCA", assay.type = "counts", ncomponents = 10)
```


```r
plotReducedDim(tse, "PCA", colour_by = "Group")
```

![(\#fig:plot-pca)PCA plot on the GlobalPatterns data set containing sample from different sources.](20_beta_diversity_files/figure-latex/plot-pca-1.png) 

As mentioned before, applicability of the different methods depends on your
sample set.

FIXME: let us switch to UMAP for the examples?


```r
tse <- runTSNE(tse, name = "TSNE", assay.type = "counts", ncomponents = 3)
```


```r
plotReducedDim(tse, "TSNE", colour_by = "Group", ncomponents = c(1:3))
```

![(\#fig:plot-tsne)t-SNE plot on the GlobalPatterns data set containing sample from different sources.](20_beta_diversity_files/figure-latex/plot-tsne-1.png) 

As a final note, `mia` provides functions for the evaluation of additional dissimilarity indices, such as:
* `calculateJSD`, `runJSD` (Jensen-Shannon divergence)
* `calculateNMDS`, `plotNMDS` (non-metric multi-dimensional scaling)
* `calculateCCA`, `runCCA` (Canonical Correspondence Analysis)
* `calculateRDA`, `runRDA` (Redundancy Analysis)
* `calculateOverlap`, `runOverlap` ()
* `calculateDPCoA`, `runDPCoA` (Double Principal Coordinate Analysis)

Redundancy analysis is similar to PCA, however, it takes into account covariates. 
It aims to maximize the variance in respect of covariates. The results shows how much
each covariate affects.


```r
# Load required packages
if(!require("vegan")){
    install.packages("vegan")
    library("vegan")
}
if(!require("stringr")){
    install.packages("stringr")
    library("stringr")
}
if(!require("knitr")){
    install.packages("knitr")
    library("knitr")
}
# Load data
data(enterotype, package="mia")
# Covariates that are being analyzed
variable_names <- c("ClinicalStatus", "Gender", "Age")

# Apply relative transform
enterotype <- transformCounts(enterotype, method = "relabundance")

# Create a formula
formula <- as.formula(paste0("assay ~ ", str_c(variable_names, collapse = " + ")) )

# # Perform RDA
rda <- calculateRDA(enterotype, assay_name = "relabundance",
                    formula = formula, distance = "bray", na.action = na.exclude)
# Get the rda object
rda <- attr(rda, "rda")
# Calculate p-value and variance for whole model
# Recommendation: use 999 permutations instead of 99
set.seed(436)
permanova <- anova.cca(rda, permutations = 99)
# Create a data.frame for results
rda_info <- as.data.frame(permanova)["Model", ]

# Calculate p-value and variance for each variable
# by = "margin" --> the order or variables does not matter
set.seed(4585)
permanova <- anova.cca(rda, by = "margin",  permutations = 99)
# Add results to data.frame
rda_info <- rbind(rda_info, permanova)

# Add info about total variance
rda_info[ , "Total variance"] <- rda_info["Model", 2] +
    rda_info["Residual", 2]

# Add info about explained variance
rda_info[ , "Explained variance"] <- rda_info[ , 2] / 
    rda_info[ , "Total variance"]

# Loop through variables, calculate homogeneity
homogeneity <- list()
# Get colDtaa
coldata <- colData(enterotype)
# Get assay
assay <- t(assay(enterotype, "relabundance"))
for( variable_name in rownames(rda_info) ){
    # If data is continuous or discrete
    if( variable_name %in% c("Model", "Residual") ||
        length(unique(coldata[[variable_name]])) /
        length(coldata[[variable_name]]) > 0.2 ){
        # Do not calculate homogeneity for continuous data
        temp <- NA
    } else{
        # Calculate homogeneity for discrete data
        # Calculate homogeneity
        set.seed(413)
        temp <- anova(
            betadisper( 
                vegdist(assay, method = "bray"),
                group = coldata[[variable_name]] ),
            permutations = permutations )["Groups", "Pr(>F)"]
    }
    # Add info to the list
    homogeneity[[variable_name]] <- temp
}
# Add homogeneity to information
rda_info[["Homogeneity p-value (NULL hyp: distinct/homogeneous --> permanova suitable)"]] <-
    homogeneity

kable(rda_info)
```


\begin{tabular}{l|r|r|r|r|r|r|l}
\hline
  & Df & SumOfSqs & F & Pr(>F) & Total variance & Explained variance & Homogeneity p-value (NULL hyp: distinct/homogeneous --> permanova suitable)\\
\hline
Model & 6 & 1.1157 & 1.940 & 0.05 & 3.991 & 0.2795 & NA\\
\hline
ClinicalStatus & 4 & 0.5837 & 1.522 & 0.15 & 3.991 & 0.1463 & 0.044277....\\
\hline
Gender & 1 & 0.1679 & 1.751 & 0.10 & 3.991 & 0.0421 & 0.522999....\\
\hline
Age & 1 & 0.5245 & 5.471 & 0.01 & 3.991 & 0.1314 & 0.000369....\\
\hline
Residual & 30 & 2.8757 & NA & NA & 3.991 & 0.7205 & NA\\
\hline
\end{tabular}


```r
# Load ggord for plotting
if(!require("ggord")){
    if(!require("devtools")){
        install.packages("devtools")
        library("devtools")
    }
    install_github("https://github.com/fawda123/ggord/")
    library("ggord")
}
if(!require("ggplot2")){
    install.packages("ggplot2")
    library("ggplot2")
}
# Since na.exclude was used, if there were rows missing information, they were 
# dropped off. Subset coldata so that it matches with rda.
coldata <- coldata[ rownames(rda$CCA$wa), ]

# Adjust names
# Get labels of vectors
vec_lab_old <- rownames(rda$CCA$biplot)

# Loop through vector labels
vec_lab <- sapply(vec_lab_old, FUN = function(name){
    # Get the variable name
    variable_name <- variable_names[ str_detect(name, variable_names) ]
    # If the vector label includes also group name
    if( !any(name %in% variable_names) ){
        # Get the group names
        group_name <- unique( coldata[[variable_name]] )[ 
        which( paste0(variable_name, unique( coldata[[variable_name]] )) == name ) ]
        # Modify vector so that group is separated from variable name
        new_name <- paste0(variable_name, " \U2012 ", group_name)
    } else{
        new_name <- name
    }
    # Add percentage how much this variable explains, and p-value
    new_name <- expr(paste(!!new_name, " (", 
                           !!format(round( rda_info[variable_name, "Explained variance"]*100, 1), nsmall = 1), 
                           "%, ",italic("P"), " = ", 
                           !!gsub("0\\.","\\.", format(round( rda_info[variable_name, "Pr(>F)"], 3), 
                                                       nsmall = 3)), ")"))

    return(new_name)
})
# Add names
names(vec_lab) <- vec_lab_old

# Create labels for axis
xlab <- paste0("RDA1 (", format(round( rda$CCA$eig[[1]]/rda$CCA$tot.chi*100, 1), nsmall = 1 ), "%)")
ylab <- paste0("RDA2 (", format(round( rda$CCA$eig[[2]]/rda$CCA$tot.chi*100, 1), nsmall = 1 ), "%)")

# Create a plot        
plot <- ggord(rda, grp_in = coldata[["ClinicalStatus"]], vec_lab = vec_lab,
              alpha = 0.5,
              size = 4, addsize = -4,
              #ext= 0.7, 
              txt = 3.5, repel = TRUE, 
              #coord_fix = FALSE
          ) + 
    # Adjust titles and labels
    guides(colour = guide_legend("ClinicalStatus"),
           fill = guide_legend("ClinicalStatus"),
           group = guide_legend("ClinicalStatus"),
           shape = guide_legend("ClinicalStatus"),
           x = guide_axis(xlab),
           y = guide_axis(ylab)) +
    theme( axis.title = element_text(size = 10) )
plot
```

![](20_beta_diversity_files/figure-latex/microbiome_RDA2-1.png)<!-- --> 

From RDA plot, we can see that only age has significant affect on microbial profile. 

## Visualizing the most dominant genus on PCoA {#pcoa-genus}

In this section we visualize most dominant genus on PCoA. A similar visualization was proposed by Salosensaari et al. [-@Salosensaari2021].


Let us agglomerate the data at a Genus level and getting the dominant taxa per sample.


```r
# Agglomerate to genus level
tse_Genus <- agglomerateByRank(tse, rank="Genus")
# Convert to relative abundances
tse_Genus <- transformCounts(tse, method = "relabundance", assay_name="counts")
# Add info on dominant genus per sample
tse_Genus <- addPerSampleDominantTaxa(tse_Genus, assay_name="relabundance", name = "dominant_taxa")
```


Performing PCoA with Bray-Curtis dissimilarity.

```r
tse_Genus <- runMDS(tse_Genus, FUN = vegan::vegdist,
              name = "PCoA_BC", assay.type = "relabundance")
```


Getting top taxa and visualizing the abundance on PCoA.


```r
# Getting the top taxa
top_taxa <- getTopTaxa(tse_Genus,top = 6, assay_name = "relabundance")

# Naming all the rest of non top-taxa as "Other"
most_abundant <- lapply(colData(tse_Genus)$dominant_taxa,
                   function(x){if (x %in% top_taxa) {x} else {"Other"}})

# Storing the previous results as a new column within colData
colData(tse_Genus)$most_abundant <- as.character(most_abundant)

# Calculating percentage of the most abundant
most_abundant_freq <- table(as.character(most_abundant))
most_abundant_percent <- round(most_abundant_freq/sum(most_abundant_freq)*100, 1)

# Retrieving the explained variance
e <- attr(reducedDim(tse_Genus, "PCoA_BC"), "eig");
var_explained <- e/sum(e[e>0])*100

# Visualization
plot <-plotReducedDim(tse_Genus,"PCoA_BC", colour_by = "most_abundant") +
  scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red"),
                      labels=paste0(names(most_abundant_percent),"(",most_abundant_percent,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       color="")
plot
```

![](20_beta_diversity_files/figure-latex/unnamed-chunk-7-1.png)<!-- --> 

Note: A 3D interactive version of the earlier plot can be found from \@ref(extras).

Similarly let's visualize and compare the sub-population.


```r
# Calculating the frequencies and percentages for both categories
freq_TRUE <- table(as.character(most_abundant[colData(tse_Genus)$Group==TRUE]))
freq_FALSE <- table(as.character(most_abundant[colData(tse_Genus)$Group==FALSE]))
percent_TRUE <- round(freq_TRUE/sum(freq_TRUE)*100, 1)
percent_FALSE <- round(freq_FALSE/sum(freq_FALSE)*100, 1)

# Visualization
plotReducedDim(tse_Genus[,colData(tse_Genus)$Group==TRUE],
                          "PCoA_BC", colour_by = "most_abundant") +
  scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red"),
                      labels=paste0(names(percent_TRUE),"(",percent_TRUE,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       title = "Group = TRUE", color="")
```

![](20_beta_diversity_files/figure-latex/unnamed-chunk-8-1.png)<!-- --> 

```r
plotReducedDim(tse_Genus[,colData(tse_Genus)$Group==FALSE],
                          "PCoA_BC", colour_by = "most_abundant") +
  scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red"),
                      labels=paste0(names(percent_FALSE),"(",percent_FALSE,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       title = "Group = FALSE", color="")
```

![](20_beta_diversity_files/figure-latex/unnamed-chunk-8-2.png)<!-- --> 



### Testing differences in community composition between sample groups

The permutational analysis of variance (PERMANOVA) [@Anderson2001] is
a widely used non-parametric multivariate method that can be used to
estimate the actual statistical significance of differences in the
observed community composition between two groups of
samples.

PERMANOVA evaluates the hypothesis that the centroids and dispersion
of the community are equivalent between the compared groups. A small
p-value indicates that the compared groups have, on average, a
different community composition.

This method is implemented in the `vegan` package in the function
[`adonis2`](https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/adonis).

**Note:**

It is recommended to `by = "margin"`. It specifies that each variable's marginal
effect is analyzed individually. 

When `by = "terms"` (the default)  the order of variables matters;
each variable is analyzed sequentially, and the result is different when more than 1 variable is
introduced and their order is differs. (Check [comparison](https://microbiome.github.io/OMA/extras.html#permanova-comparison))

We can perform PERMANOVA with `adonis2` function or by first performing distance-based 
redundancy analysis (dbRDA), and then applying permutational test for result of 
redundancy analysis. Advantage of the latter approach is that by doing so we can get 
coefficients: how much each taxa affect to the result.


```r
if( !require(vegan) ){
    BiocManager::install("vegan")
    library("vegan")
}
# Agglomerate data to Species level
tse <- agglomerateByRank(tse, rank = "Species")

# Set seed for reproducibility
set.seed(1576)
# We choose 99 random permutations. Consider applying more (999 or 9999) in your
# analysis. 
permanova <- adonis2(t(assay(tse,"relabundance")) ~ Group,
                     by = "margin", # each term (here only 'Group') analyzed individually
                     data = colData(tse),
                     method = "euclidean",
                     permutations = 99)

# Set seed for reproducibility
set.seed(1576)
# Perform dbRDA
dbrda <- dbrda(t(assay(tse,"relabundance")) ~ Group, 
               data = colData(tse))
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term (here only 'Group') analyzed individually
                        method = "euclidean",
                        permutations = 99)

# Get p-values
p_values <- c( permanova["Group", "Pr(>F)"], permanova2["Group", "Pr(>F)"] )
p_values <-as.data.frame(p_values)
rownames(p_values) <- c("adonis2", "dbRDA+anova.cca")
p_values
```

```
##                 p_values
## adonis2             0.02
## dbRDA+anova.cca     0.02
```

As we can see, the community composition is significantly different
between the groups (p < 0.05), and these two methods give equal p-values.

Let us visualize the model coefficients for species that exhibit the
largest differences between the groups. This gives some insights into
how the groups tend to differ from each other in terms of community
composition.


```r
# Add taxa info
sppscores(dbrda) <- t(assay(tse,"relabundance"))
# Get coefficients
coef <- dbrda$CCA$v
# Get the taxa with biggest weights
top.coef <- head( coef[rev(order(abs(coef))), , drop = FALSE], 20)
# Sort weights in increasing order
top.coef <- top.coef[ order(top.coef), ]
# Get top names
top_names <- names(top.coef)[ order(abs(top.coef), decreasing = TRUE) ]
```




```r
ggplot(data.frame(x = top.coef,
                  y = factor(names(top.coef),
                                      unique(names(top.coef)))),
        aes(x = x, y = y)) +
    geom_bar(stat="identity") +
    labs(x="",y="",title="Top Taxa") +
    theme_bw()
```

![](20_beta_diversity_files/figure-latex/plot-top-coef-anova-1.png)<!-- --> 

In the above example, the largest differences between the two groups
can be attributed to _Genus:Bacteroides_ (elevated in the first
group) and _Family:Ruminococcaceae_ (elevated in the second
group), and many other co-varying species.



### Checking the homogeneity condition 

It is important to note that the application of PERMANOVA assumes
homogeneous group dispersions (variances). This can be tested with the
PERMDISP2 method [@Anderson2006] by using the same assay and distance
method than in PERMANOVA.


```r
anova( betadisper(vegdist(t(assay(tse, "counts"))), colData(tse)$Group) )
```

```
## Analysis of Variance Table
## 
## Response: Distances
##           Df Sum Sq Mean Sq F value  Pr(>F)    
## Groups     1 0.2385  0.2385     103 3.6e-10 ***
## Residuals 24 0.0554  0.0023                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

If the groups have similar dispersion, PERMANOVA can be seen as an
appropriate choice for comparing community compositions.


## Further reading

 - [How to extract information from clusters](http://bioconductor.org/books/release/OSCA/clustering.html)

 - Chapter \@ref(clustering) on community typing

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
 [1] ggord_1.1.7                    knitr_1.42                    
 [3] stringr_1.5.0                  patchwork_1.1.2               
 [5] scater_1.28.0                  scuttle_1.10.0                
 [7] ggplot2_3.4.2                  vegan_2.6-4                   
 [9] lattice_0.21-8                 permute_0.9-7                 
[11] mia_1.9.2                      MultiAssayExperiment_1.26.0   
[13] TreeSummarizedExperiment_2.1.4 Biostrings_2.68.0             
[15] XVector_0.40.0                 SingleCellExperiment_1.22.0   
[17] SummarizedExperiment_1.30.0    Biobase_2.60.0                
[19] GenomicRanges_1.52.0           GenomeInfoDb_1.36.0           
[21] IRanges_2.34.0                 S4Vectors_0.38.0              
[23] BiocGenerics_0.46.0            MatrixGenerics_1.12.0         
[25] matrixStats_0.63.0-9003        BiocStyle_2.28.0              
[27] rebook_1.9.0                  

loaded via a namespace (and not attached):
 [1] DBI_1.1.3                   bitops_1.0-7               
 [3] gridExtra_2.3               CodeDepends_0.6.5          
 [5] rlang_1.1.1                 magrittr_2.0.3             
 [7] compiler_4.3.0              RSQLite_2.3.1              
 [9] mgcv_1.8-42                 dir.expiry_1.8.0           
[11] DelayedMatrixStats_1.22.0   vctrs_0.6.2                
[13] reshape2_1.4.4              pkgconfig_2.0.3            
[15] crayon_1.5.2                fastmap_1.1.1              
[17] labeling_0.4.2              utf8_1.2.3                 
[19] rmarkdown_2.21              graph_1.78.0               
[21] ggbeeswarm_0.7.2            DirichletMultinomial_1.42.0
[23] purrr_1.0.1                 bit_4.0.5                  
[25] xfun_0.39                   zlibbioc_1.46.0            
[27] cachem_1.0.7                beachmat_2.16.0            
[29] jsonlite_1.8.4              blob_1.2.4                 
[31] highr_0.10                  DelayedArray_0.25.0        
[33] BiocParallel_1.34.0         cluster_2.1.4              
[35] irlba_2.3.5.1               parallel_4.3.0             
[37] R6_2.5.1                    stringi_1.7.12             
[39] Rcpp_1.0.10                 bookdown_0.33              
[41] DECIPHER_2.28.0             splines_4.3.0              
[43] Matrix_1.5-4                tidyselect_1.2.0           
[45] yaml_2.3.7                  viridis_0.6.2              
[47] codetools_0.2-19            tibble_3.2.1               
[49] plyr_1.8.8                  withr_2.5.0                
[51] treeio_1.24.0               Rtsne_0.16                 
[53] evaluate_0.20               pillar_1.9.0               
[55] BiocManager_1.30.20         filelock_1.0.2             
[57] generics_0.1.3              RCurl_1.98-1.12            
[59] sparseMatrixStats_1.12.0    munsell_0.5.0              
[61] scales_1.2.1                tidytree_0.4.2             
[63] glue_1.6.2                  lazyeval_0.2.2             
[65] tools_4.3.0                 BiocNeighbors_1.18.0       
[67] ScaledMatrix_1.7.1          XML_3.99-0.14              
[69] cowplot_1.1.1               grid_4.3.0                 
[71] tidyr_1.3.0                 ape_5.7-1                  
[73] colorspace_2.1-0            nlme_3.1-162               
[75] GenomeInfoDbData_1.2.10     beeswarm_0.4.0             
[77] BiocSingular_1.16.0         vipor_0.4.5                
[79] cli_3.6.1                   rsvd_1.0.5                 
[81] fansi_1.0.4                 viridisLite_0.4.1          
[83] dplyr_1.1.2                 gtable_0.3.3               
[85] yulab.utils_0.0.6           digest_0.6.31              
[87] ggrepel_0.9.3               farver_2.1.1               
[89] decontam_1.20.0             memoise_2.0.1              
[91] htmltools_0.5.5             lifecycle_1.0.3            
[93] bit64_4.0.5                 MASS_7.3-59                
```
</div>
