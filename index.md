--- 
title: "Orchestrating Microbiome Analysis"
documentclass: book
site: bookdown::bookdown_site
bibliography: [book.bib]
biblio-style: apalike
link-citations: yes
github-repo: microbiome/OMA
always_allow_html: yes
---



---
date: "**Authors:** Leo Lahti [aut], Sudarshan Shetty [aut], Felix GM Ernst [aut, cre]<br/>
  **Version:** 0.98.9<br/>
  **Modified:** 2021-09-10<br/>
  **Compiled:** 2022-02-27<br/>
  **Environment:** R version 4.1.2 (2021-11-01), Bioconductor 3.14<br/>
  **License:** CC BY-NC-SA 3.0 US<br/>
  **Copyright:** <br/>
  **Source:** https://github.com/microbiome/OMA"
url: "https://github.com/microbiome/OMA"
---





# Welcome {-}


<a href="https://bioconductor.org"><img src="https://github.com/Bioconductor/BiocStickers/raw/master/Bioconductor/Bioconductor-serial.gif" width="200" alt="Bioconductor Sticker" align="right" style="margin: 0 1em 0 1em" /></a>


You are reading the online book, [**Orchestrating Microbiome Analysis
with R and Bioconductor**](microbiome.github.io/OMA) [@OMA], where we
walk through common strategies and workflows in microbiome data
science.

The book shows through concrete examples how you can take advantage of
the latest developments in R/Bioconductor for the manipulation,
analysis, and reproducible reporting of hierarchical and heterogeneous
microbiome profiling data sets. The book was borne out of necessity,
while updating microbiome analysis tools to work with Bioconductor
classes that provide support for multi-modal data collections. Many of
these techniques are generic and widely applicable in other contexts
as well.

This work has been heavily influenced by other similar resources, in
particular the Orchestrating Single-Cell Analysis with Bioconductor
[@Amezquita2020natmeth], [phyloseq
tutorials](http://joey711.github.io/phyloseq/tutorials-index)
[@Callahan2016] and [microbiome
tutorials](https://microbiome.github.io/tutorials/) [@Shetty2019].
This book extends these resources to teach the grammar of Bioconductor
workflows in the context of microbiome data science.  As such, it
supports the adoption of general skills in the analysis of large,
hierarchical, and multi-modal data collections. We focus on microbiome
analysis tools, including entirely new, partially updated as well as
previously established methods.

This online resource and its associated ecosystem of microbiome data
science tools are a result of a community-driven development process,
and welcoming new contributors. Several individuals have
[contributed](https://github.com/microbiome/OMA/graphs/contributors)
methods, workflows and improvements as acknowledged in the
Introduction. You can find more information on how to find us online
and join the developer community through the project homepage at
[microbiome.github.io](https://microbiome.github.io). This online
resource has been written in RMarkdown with the bookdown R
package. The material is **free to use** with the [Creative Commons
Attribution-NonCommercial
3.0](https://creativecommons.org/licenses/by-nc/3.0/us/) License.


--------------



