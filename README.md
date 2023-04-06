# Orchestrating Microbiome Analysis Book <img src="general/figures/mia_logo.png" align="right" width="120" />

## Overview (second version)

This is a reference cookbook for performing **Microbiome Analysis** with 
Bioconductor in R. This is a book based on R Markdown and **bookdown** 
(https://github.com/rstudio/bookdown).

## Deployment

For now the book is deployed to GitHub Pages from GitHub Actions. 
The book can also be built locally:

```
bookdown::render_book("index.Rmd", "bookdown::gitbook")
```

You can view the book also with:

```
bookdown::serve_book()
```


To install the necessary dependencies to build to book, run

```
BiocManager::install(remotes::local_package_deps(dependencies=TRUE))
```

## Development and Contribtions

To contribute reports, follow the Git flow procedure:

1. Fork the project
2. Create a new branch
3. Commit changes to the new branch
4. Create a pull request (PR) back to the original repo
5. Fix and discuss issues in the review process

Please note that chapters should be independent of each other.

### Adding new sections

- Create the relevant Rmd file; follow the numbering logic
- Add it also to the list in file [_bookdown.yml](_bookdown.yml). 


# Code of conduct

Please note that the OMA project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
