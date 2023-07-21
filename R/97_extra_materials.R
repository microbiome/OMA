## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ----permanova_import, warning = FALSE, message = FALSE-----------------------
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


## ----permanova_prep, message = FALSE, warning = FALSE-------------------------
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


## ----permanova_loop, message = FALSE, warning = FALSE-------------------------
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


## ----permanova_table, message = FALSE, warning = FALSE------------------------

df <- terms_df %>%
  dplyr::inner_join(margin_df, by = "Formula", suffix = c(" (terms)", " (margin)"))

knitr::kable(df)


## ---- message=FALSE, warning=FALSE--------------------------------------------
if (!require(fido)){
  # installing the fido package
  devtools::install_github("jsilve24/fido")
}


## ---- message=FALSE, warning=FALSE--------------------------------------------
library(fido)


## ---- message=FALSE, warning=FALSE, eval=FALSE--------------------------------
## library(microbiomeDataSets)
## tse <- SprockettTHData()


## ---- message=FALSE, warning=FALSE, echo=FALSE--------------------------------
# saveRDS(tse, file="data/SprockettTHData.Rds")
# Hidden reading of the saved data
tse <- readRDS("data/SprockettTHData.Rds")


## ---- message=FALSE, warning=FALSE--------------------------------------------
library(mia)
cov_names <- c("Sex","Age_Years","Delivery_Mode")
na_counts <- apply(is.na(colData(tse)[,cov_names]), 2, sum)
na_summary<-as.data.frame(na_counts,row.names=cov_names)


## ---- message=FALSE, warning=FALSE--------------------------------------------
tse <- tse[ , !is.na(colData(tse)$Delivery_Mode) ]
tse <- tse[ , !is.na(colData(tse)$Age_Years) ]


## ---- message=FALSE, warning=FALSE--------------------------------------------
tse_phylum <- agglomerateByRank(tse, "Phylum")


## ---- message=FALSE, warning=FALSE--------------------------------------------
Y <- assays(tse_phylum)$counts
# design matrix
# taking 3 covariates
sample_data<-as.data.frame(colData(tse_phylum)[,cov_names])
X <- t(model.matrix(~Sex+Age_Years+Delivery_Mode,data=sample_data))


## ---- message=FALSE, warning=FALSE--------------------------------------------
n_taxa<-nrow(Y)
upsilon <- n_taxa+3
Omega <- diag(n_taxa)
G <- cbind(diag(n_taxa-1), -1)
Xi <- (upsilon-n_taxa)*G%*%Omega%*%t(G)
Theta <- matrix(0, n_taxa-1, nrow(X))
Gamma <- diag(nrow(X))


## ---- message=FALSE, warning=FALSE--------------------------------------------
priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)
names_covariates(priors) <- rownames(X)
plot(priors, pars="Lambda") + ggplot2::xlim(c(-5, 5))


## ---- message=FALSE, warning=FALSE--------------------------------------------
priors$Y <- Y 
posterior <- refit(priors, optim_method="adam", multDirichletBoot=0.5) #calcGradHess=FALSE


## ---- message=FALSE, warning=FALSE--------------------------------------------
ppc_summary(posterior)


## ---- message=FALSE, warning=FALSE--------------------------------------------
names_categories(posterior) <- rownames(Y)
plot(posterior,par="Lambda",focus.cov=rownames(X)[2:4])


## ---- message=FALSE, warning=FALSE--------------------------------------------
plot(posterior, par="Lambda", focus.cov = rownames(X)[c(2,4)])


## ---- message=FALSE, warning=FALSE--------------------------------------------
# Installing required packages
if (!require(rgl)){
  BiocManager::install("rgl")  
}
if (!require(plotly)){
  BiocManager::install("plotly")  
}


## ----setup2, warning=FALSE, message=FALSE-------------------------------------
library(knitr)
library(rgl)
knitr::knit_hooks$set(webgl = hook_webgl)


## ---- message=FALSE, warning=FALSE--------------------------------------------
# Installing the package
if (!require(curatedMetagenomicData)){
  BiocManager::install("curatedMetagenomicData")  
}
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

tse_Genus <- agglomerateByRank(tse, rank="genus")
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

