## ---- echo=FALSE---------------------------------------
.gh_url <- file.path('https://github.com', rmarkdown::metadata[['github-repo']])


## ----include=FALSE-------------------------------------
library(Cairo)

# global knitr options
knitr::opts_chunk$set(
  fig.width=10,
  dpi=300,
  dev = "png",
  dev.args = list(type = "cairo-png")
)


## // This block adds image to the front page

## title=document.getElementById('header');

## title.innerHTML = title.innerHTML +

## 

## '<img src="https://user-images.githubusercontent.com/60338854/128359392\

## -6feef8df-30e9-4ea0-ae3b-4bb619d746ed.png" alt="Microbiome" width="50%"/>' +

## 

## '<p style="font-size:12px">Figure source: Moreno-Indias <i>et al</i>. (2021) \

## <a href="https://doi.org/10.3389/fmicb.2021.635781">Statistical and \

## Machine Learning Techniques in Human Microbiome Studies: Contemporary \

## Challenges and Solutions</a>. Frontiers in Microbiology 12:11.</p>'

