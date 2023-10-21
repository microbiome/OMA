## ----setup-example, eval = FALSE---------------------------------------------------------------------------------------
## # setup
## library(mia)
## data("GlobalPatterns", package = "mia")
## tse <- GlobalPatterns
## 
## # this line sets some options for all the chunks (global chunk options)
## knitr::opts_chunk$set(message = FALSE, warning = FALSE)


## ----fig-example, eval = FALSE-----------------------------------------------------------------------------------------
## # fig-box
## boxplot(colSums(assay(tse)) ~ tse$SampleType)


## ----tbl-example, eval = FALSE-----------------------------------------------------------------------------------------
## # tbl-coldata
## knitr::kable(head(colData(tse)))

