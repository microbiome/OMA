# Extra material {#extras}

## PERMANOVA comparison

Here we present two possible uses of the `adonis2` function which performs PERMANOVA. The
optional argument `by` has an effect on the statistical outcome, so its two options are
compared here.


```r
# import necessary packages
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

Let us load the _enterotype_ TSE object and run PERMANOVA for
different orders of three variables with two different approaches:
`by = "margin"` or `by = "terms"`.



```r
# load and prepare data
library(mia)
data("enterotype", package="mia")
enterotype <- transformCounts(enterotype, method = "relabundance")
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


