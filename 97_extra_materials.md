# Extra material {#extras}

## Interactive 3D Plots


```r
# Installing required packages
if (!require(rgl)){
  BiocManager::install("rgl")  
}
if (!require(plotly)){
  BiocManager::install("plotly")  
}
```


```r
library(knitr)
library(rgl)
knitr::knit_hooks$set(webgl = hook_webgl)
```


In this section we make a 3D version of the earlier [ Visualizing the most dominant genus on PCoA](https://microbiome.github.io/OMA/microbiome-exploration.html#visualizing-the-most-dominant-genus-on-pcoa), with the help of the [plotly](https://plotly.com/r/) R package.


```r
# Importing necessary libraries
library(curatedMetagenomicData)
library(dplyr)
library(DT)
library(mia)
library(scater)

# Querying the data
tse <- sampleMetadata |>
    filter(age >= 18) |> # taking only data of age 18 or above
    filter(!is.na(alcohol)) |> # excluding missing values
    select(where(~ !all(is.na(.x)))) |>
    returnSamples("relative_abundance")

tse_Genus <- agglomerateByRank(tse, rank="genus")
tse_Genus <- addPerSampleDominantTaxa(tse_Genus, assay_name="relative_abundance", name = "dominant_taxa")

# Performing PCoA with Bray-Curtis dissimilarity.
tse_Genus <- runMDS(tse_Genus, FUN = vegan::vegdist, ncomponents = 3,
              name = "PCoA_BC", exprs_values = "relative_abundance")

# Getting the 6 top taxa
top_taxa <- getTopTaxa(tse_Genus,top = 6, assay_name = "relative_abundance")

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
```

Interactive 3D visualization of the most dominant genus on PCoA.

Note that labels at legend can be used to visualize one or more Genus
separately (double click to isolate one from the others, or toggle to
select multiple ones).


```r
library(plotly)

# 3D Visualization
reduced_data  <- as.data.frame(reducedDim(tse_Genus)[,])
names(reduced_data) <- c("PC1","PC2","PC3")
plot_ly(reduced_data, x=~PC1,y=~PC2,z=~PC3) %>%
  add_markers(color=sapply(strsplit(colData(tse_Genus)$most_abundant, "_"), tail, 1), size=5,
              colors=c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red")) %>%
  layout(scene=list(xaxis=list(title = paste("PC1 (",round(var_explained[1],1),"%)")),
                    yaxis=list(title = paste("PC2 (",round(var_explained[2],1),"%)")),
                    zaxis=list(title = paste("PC3 (",round(var_explained[3],1),"%)"))))
```

![](97_extra_materials_files/figure-latex/test-rgl-1.pdf)<!-- --> 


## PERMANOVA comparison


```r
library("vegan")
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.6-4
```

```r
data("enterotype")
enterotype <- transformSamples(enterotype, method = "relabundance")
```

```
## Warning: useNames = NA is deprecated. Instead, specify either useNames = TRUE
## or useNames = FALSE.
```

```r
# Drop those samples that do not have meta dtaa
enterotype <- 
    enterotype[ , !rowSums(is.na(colData(enterotype)[, c("Nationality", "Gender","ClinicalStatus")]) > 0 ) ]

# Multiple variables, by = "margin"
set.seed(75)
adonis2(t(assay(enterotype ,"relabundance")) ~ Nationality + Gender + ClinicalStatus,
        by = "margin",
        data = colData(enterotype),
        permutations = 99)
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 99
## 
## adonis2(formula = t(assay(enterotype, "relabundance")) ~ Nationality + Gender + ClinicalStatus, data = colData(enterotype), permutations = 99, by = "margin")
##                Df SumOfSqs      R2      F Pr(>F)  
## Nationality     4   0.7933 0.19338 2.0894   0.04 *
## Gender          1   0.1230 0.02999 1.2963   0.29  
## ClinicalStatus  3   0.2412 0.05879 0.8469   0.53  
## Residual       29   2.7527 0.67101                
## Total          38   4.1023 1.00000                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Multiple variables, by = "margin"
set.seed(75)
adonis2(t(assay(enterotype ,"relabundance")) ~ ClinicalStatus + Nationality + Gender,
        by = "margin",
        data = colData(enterotype ),
        permutations = 99)
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 99
## 
## adonis2(formula = t(assay(enterotype, "relabundance")) ~ ClinicalStatus + Nationality + Gender, data = colData(enterotype), permutations = 99, by = "margin")
##                Df SumOfSqs      R2      F Pr(>F)  
## ClinicalStatus  3   0.2412 0.05879 0.8469   0.53  
## Nationality     4   0.7933 0.19338 2.0894   0.04 *
## Gender          1   0.1230 0.02999 1.2963   0.29  
## Residual       29   2.7527 0.67101                
## Total          38   4.1023 1.00000                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Multiple variables, default: by = "terms"
set.seed(75)
adonis2(t(assay(enterotype ,"relabundance"))  ~ Nationality + Gender + ClinicalStatus,
        data = colData(enterotype ),
        permutations = 99)
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 99
## 
## adonis2(formula = t(assay(enterotype, "relabundance")) ~ Nationality + Gender + ClinicalStatus, data = colData(enterotype), permutations = 99)
##                Df SumOfSqs      R2      F Pr(>F)  
## Nationality     5   1.0046 0.24490 2.1168   0.04 *
## Gender          1   0.1038 0.02530 1.0934   0.39  
## ClinicalStatus  3   0.2412 0.05879 0.8469   0.53  
## Residual       29   2.7527 0.67101                
## Total          38   4.1023 1.00000                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Multiple variables, default: by = "terms"
set.seed(75)
adonis2(t(assay(enterotype ,"relabundance"))  ~ ClinicalStatus + Nationality + Gender,
        data = colData(enterotype),
        permutations = 99)
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 99
## 
## adonis2(formula = t(assay(enterotype, "relabundance")) ~ ClinicalStatus + Nationality + Gender, data = colData(enterotype), permutations = 99)
##                Df SumOfSqs      R2      F Pr(>F)  
## ClinicalStatus  4   0.5000 0.12189 1.3169   0.20  
## Nationality     4   0.7265 0.17710 1.9135   0.05 *
## Gender          1   0.1230 0.02999 1.2963   0.29  
## Residual       29   2.7527 0.67101                
## Total          38   4.1023 1.00000                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Bayesian Multinomial Logistic-Normal Models

Analysis using such model could be performed with the function `pibble` from the `fido`
package, wihch is in form of a Multinomial Logistic-Normal Linear Regression model; see
[vignette](https://jsilve24.github.io/fido/articles/introduction-to-fido.html) of package.


The following presents such an exemplary analysis based on the 
data of @Sprockett2020 available
through `microbiomeDataSets` package.



```r
if (!require(fido)){
  # installing the fido package
  devtools::install_github("jsilve24/fido")
}
```

Loading the libraries and importing data:


```r
library(fido)
```


```r
library(microbiomeDataSets)
tse <- SprockettTHData()
```




We pick three covariates ("Sex","Age_Years","Delivery_Mode") during this
analysis as an example, and beforehand we check for missing data:



```r
cov_names <- c("Sex","Age_Years","Delivery_Mode")
na_counts <- apply(is.na(colData(tse)[,cov_names]), 2, sum)
na_summary<-as.data.frame(na_counts,row.names=cov_names)
```

We drop missing values of the covariates:


```r
tse <- tse[ , !is.na(colData(tse)$Delivery_Mode) ]
tse <- tse[ , !is.na(colData(tse)$Age_Years) ]
```

We agglomerate microbiome data to Phylum:


```r
tse_phylum <- agglomerateByRank(tse, "Phylum")
```

We extract the counts assay and covariate data to build the model
matrix:


```r
Y <- assays(tse_phylum)$counts
# design matrix
# taking 3 covariates
sample_data<-as.data.frame(colData(tse_phylum)[,cov_names])
X <- t(model.matrix(~Sex+Age_Years+Delivery_Mode,data=sample_data))
```

Building the parameters for the `pibble` call to build the model; see more at [vignette](https://jsilve24.github.io/fido/articles/introduction-to-fido.html):


```r
n_taxa<-nrow(Y)
upsilon <- n_taxa+3
Omega <- diag(n_taxa)
G <- cbind(diag(n_taxa-1), -1)
Xi <- (upsilon-n_taxa)*G%*%Omega%*%t(G)
Theta <- matrix(0, n_taxa-1, nrow(X))
Gamma <- diag(nrow(X))
```

Automatically initializing the priors and visualizing their distributions:


```r
priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)
names_covariates(priors) <- rownames(X)
plot(priors, pars="Lambda") + ggplot2::xlim(c(-5, 5))
```

![](97_extra_materials_files/figure-latex/unnamed-chunk-13-1.pdf)<!-- --> 

Estimating the posterior by including our response data `Y`.
Note: Some computational failures could occur (see [discussion](https://github-wiki-see.page/m/jsilve24/fido/wiki/Frequently-Asked-Questions))
the arguments `multDirichletBoot` `calcGradHess` could be passed in such case.


```r
priors$Y <- Y 
posterior <- refit(priors, optim_method="adam", multDirichletBoot=0.5) #calcGradHess=FALSE
```

Printing a summary about the posterior:


```r
ppc_summary(posterior)
```

```
## Proportions of Observations within 95% Credible Interval: 0.9978261
```
Plotting the summary of the posterior distributions of the regression parameters:


```r
names_categories(posterior) <- rownames(Y)
plot(posterior,par="Lambda",focus.cov=rownames(X)[2:4])
```

![](97_extra_materials_files/figure-latex/unnamed-chunk-16-1.pdf)<!-- --> 

Taking a closer look at "Sex" and "Delivery_Mode":


```r
plot(posterior, par="Lambda", focus.cov = rownames(X)[c(2,4)])
```

![](97_extra_materials_files/figure-latex/unnamed-chunk-17-1.pdf)<!-- --> 
