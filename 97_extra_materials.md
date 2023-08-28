# (PART) Appendix {-}

# Extra material {#extras}

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


## PERMANOVA comparison {#compare-permanova}

Here we present two possible uses of the `adonis2` function which performs PERMANOVA. The
optional argument `by` has an effect on the statistical outcome, so its two options are
compared here.


```r
# import necessary packages
library(gtools)
library(purrr)
library(vegan)
library(gtools)
library(purrr)
```

Let us load the _enterotype_ TSE object and run PERMANOVA for
different orders of three variables with two different approaches:
`by = "margin"` or `by = "terms"`.



```r
# load and prepare data
library(mia)
data("enterotype", package="mia")
enterotype <- transformAssay(enterotype, method = "relabundance")
# drop samples missing meta data
enterotype <- enterotype[ , !rowSums(is.na(colData(enterotype)[, c("Nationality", "Gender", "ClinicalStatus")]) > 0)]
# define variables and list all possible combinations
vars <- c("Nationality", "Gender", "ClinicalStatus")
var_perm <- permutations(n = 3, r = 3, vars)
formulas <- apply(var_perm, 1, function(row) purrr::reduce(row, function(x, y) paste(x, "+", y)))
# create empty data.frames for further storing p-values
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
  
  # generate temporary formula (i.e. "assay ~ ClinicalStatus + Nationality + Gender")
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

  # extract p-values
  terms_p <- with_terms[["Pr(>F)"]]
  terms_p <- terms_p[!is.na(terms_p)]
  margin_p <- with_margin[["Pr(>F)"]]
  margin_p <- margin_p[!is.na(margin_p)]
  
  # store p-values into data.frames
  for (col_idx in 1:ncol(var_perm)) {
    
    terms_df[var_perm[row_idx, col_idx]][row_idx, ] <- terms_p[col_idx]
    margin_df[var_perm[row_idx, col_idx]][row_idx, ] <- margin_p[col_idx]
    
  }
  
}
```




The following table displays the p-values for the three variables
ClinicalStatus, Gender and Nationality obtained by PERMANOVA with
`adonis2`. Note that the p-values remain identical when `by =
"margin"`, but change with the order of the variables in the
formula when `by = "terms"` (default).



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

Analysis using such model could be performed with the function
`pibble` from the `fido` package, wihch is in form of a Multinomial
Logistic-Normal Linear Regression model; see
[vignette](https://jsilve24.github.io/fido/articles/introduction-to-fido.html)
of package.


The following presents such an exemplary analysis based on the 
data of @Sprockett2020 available
through `microbiomeDataSets` package.



```r
library(fido)
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
library(mia)
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
tse_phylum <- mergeFeaturesByRank(tse, "Phylum")
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

![](97_extra_materials_files/figure-latex/unnamed-chunk-10-1.pdf)<!-- --> 

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
## Proportions of Observations within 95% Credible Interval: 0.9979
```
Plotting the summary of the posterior distributions of the regression parameters:


```r
names_categories(posterior) <- rownames(Y)
plot(posterior,par="Lambda",focus.cov=rownames(X)[2:4])
```

![](97_extra_materials_files/figure-latex/unnamed-chunk-13-1.pdf)<!-- --> 

Taking a closer look at "Sex" and "Delivery_Mode":


```r
plot(posterior, par="Lambda", focus.cov = rownames(X)[c(2,4)])
```

![](97_extra_materials_files/figure-latex/unnamed-chunk-14-1.pdf)<!-- --> 


## Interactive 3D Plots


```r
# Installing libraryd packages
library(rgl)
library(plotly)
```


```r
library(knitr)
library(rgl)
knitr::knit_hooks$set(webgl = hook_webgl)
```


In this section we make a 3D version of the earlier  Visualizing the most dominant genus on PCoA (see \@ref(quality-control)), with the help of the plotly [@Sievert2020].


```r
# Installing the package
library(curatedMetagenomicData)
# Importing necessary libraries
library(curatedMetagenomicData)
library(dplyr)
library(DT)
library(mia)
library(scater)

# Querying the data
tse <- sampleMetadata %>%
    filter(age >= 18) %>% # taking only data of age 18 or above
    filter(!is.na(alcohol)) %>% # excluding missing values
    returnSamples("relative_abundance")

tse_Genus <- mergeFeaturesByRank(tse, rank="genus")
tse_Genus <- addPerSampleDominantFeatures(tse_Genus,assay.type="relative_abundance", name = "dominant_taxa")

# Performing PCoA with Bray-Curtis dissimilarity.
tse_Genus <- runMDS(tse_Genus, FUN = vegan::vegdist, ncomponents = 3,
              name = "PCoA_BC", assay.type = "relative_abundance")

# Getting the 6 top taxa
top_taxa <- getTopFeatures(tse_Genus,top = 6, assay.type = "relative_abundance")

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

