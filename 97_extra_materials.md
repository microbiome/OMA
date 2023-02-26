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
if (!require(gtools)){
  install.packages("gtools")  
}
if (!require(purrr)){
  install.packages("purrr")  
}

library(vegan)
library(gtools)
library(purrr)
```

Here we present two possible uses of the `adonis2` function which performs PERMANOVA. The
optional argument `by` has an effect on the statistical outcome, so the two possibilities are
compared in this section.


```r
data("enterotype")
enterotype <- transformCounts(enterotype, method = "relabundance")

# drop samples missing meta data
enterotype <- enterotype[ , !rowSums(is.na(colData(enterotype)[ , c("Nationality", "Gender","ClinicalStatus")]) > 0 )]

vars <- c("Nationality", "Gender", "ClinicalStatus")
by_args <- factor(c("terms", "margin"))

var_perm <- permutations(n = 3, r = 3, vars)
formulas <- apply(var_perm, 1, function(row) purrr::reduce(row, function(x, y) paste(x, "+", y)))

terms_df <- data.frame("Formula" = formulas,
                       "ClinicalStatus" = rep(0, 6),
                       "Gender" = rep(0, 6),
                       "Nationality" = rep(0, 6))

margin_df <- data.frame("Formula" = formulas,
                        "ClinicalStatus" = rep(0, 6),
                        "Gender" = rep(0, 6),
                        "Nationality" = rep(0, 6))
```



```r
for (row_idx in 1:nrow(var_perm)) {
  
  tmp_formula <- purrr::reduce(var_perm[row_idx, ], function(x, y) paste(x, "+", y))
  tmp_formula <- as.formula(paste0('t(assay(enterotype, "relabundance")) ~ ',
                            tmp_formula))

  # multiple variables, default: by = "terms"
  set.seed(75)
  with_terms <- adonis2(tmp_formula, 
                by = "terms",
                data = colData(enterotype),
                permutations = 99)
  
  # multiple variables, by = "margin"
  set.seed(75)
  with_margin <- adonis2(tmp_formula, 
                 by = "margin",
                 data = colData(enterotype),
                 permutations = 99)

  terms_p <- with_terms[["Pr(>F)"]]
  terms_p <- terms_p[!is.na(terms_p)]
  
  margin_p <- with_margin[["Pr(>F)"]]
  margin_p <- margin_p[!is.na(margin_p)]
  
  for (col_idx in 1:ncol(var_perm)) {
    
    terms_df[var_perm[row_idx, col_idx]][row_idx, ] <- terms_p[col_idx]
    margin_df[var_perm[row_idx, col_idx]][row_idx, ] <- margin_p[col_idx]
    
  }
  
}
```

The following table displays the p-values for the three variables ClinicalStatus, Gender and
Nationality obtained by PERMANOVA with `adonis2`. As you can see, p-values remain identical when
`by = "margin"`, but they change with the order of the variables in the formula when `by = terms` (default).


```r
df <- terms_df %>%
  dplyr::inner_join(margin_df, by = "Formula", suffix = c(" (terms)", " (margin)"))

knitr::kable(df)
```


\begin{tabular}{l|r|r|r|r|r|r}
\hline
Formula & ClinicalStatus (terms) & Gender (terms) & Nationality (terms) & ClinicalStatus (margin) & Gender (margin) & Nationality (margin)\\
\hline
ClinicalStatus + Gender + Nationality & 0.20 & 0.70 & 0.04 & 0.53 & 0.29 & 0.04\\
\hline
ClinicalStatus + Nationality + Gender & 0.20 & 0.29 & 0.05 & 0.53 & 0.29 & 0.04\\
\hline
Gender + ClinicalStatus + Nationality & 0.17 & 0.79 & 0.04 & 0.53 & 0.29 & 0.04\\
\hline
Gender + Nationality + ClinicalStatus & 0.53 & 0.79 & 0.03 & 0.53 & 0.29 & 0.04\\
\hline
Nationality + ClinicalStatus + Gender & 0.61 & 0.29 & 0.04 & 0.53 & 0.29 & 0.04\\
\hline
Nationality + Gender + ClinicalStatus & 0.53 & 0.39 & 0.04 & 0.53 & 0.29 & 0.04\\
\hline
\end{tabular}

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

![](97_extra_materials_files/figure-latex/unnamed-chunk-11-1.pdf)<!-- --> 

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
## Proportions of Observations within 95% Credible Interval: 0.9979296
```
Plotting the summary of the posterior distributions of the regression parameters:


```r
names_categories(posterior) <- rownames(Y)
plot(posterior,par="Lambda",focus.cov=rownames(X)[2:4])
```

![](97_extra_materials_files/figure-latex/unnamed-chunk-14-1.pdf)<!-- --> 

Taking a closer look at "Sex" and "Delivery_Mode":


```r
plot(posterior, par="Lambda", focus.cov = rownames(X)[c(2,4)])
```

![](97_extra_materials_files/figure-latex/unnamed-chunk-15-1.pdf)<!-- --> 
