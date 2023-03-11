# Extra material {#extra_bayes}

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

![](97_extra_materials_bayes_files/figure-latex/unnamed-chunk-10-1.pdf)<!-- --> 

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
## Proportions of Observations within 95% Credible Interval: 0.9938923
```
Plotting the summary of the posterior distributions of the regression parameters:


```r
names_categories(posterior) <- rownames(Y)
plot(posterior,par="Lambda",focus.cov=rownames(X)[2:4])
```

![](97_extra_materials_bayes_files/figure-latex/unnamed-chunk-13-1.pdf)<!-- --> 

Taking a closer look at "Sex" and "Delivery_Mode":


```r
plot(posterior, par="Lambda", focus.cov = rownames(X)[c(2,4)])
```

![](97_extra_materials_bayes_files/figure-latex/unnamed-chunk-14-1.pdf)<!-- --> 
