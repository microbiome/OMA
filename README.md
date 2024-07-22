<!-- badges: start -->
📦 [Repo](https://github.com/js2264/OMA) [![rworkflows](https://img.shields.io/github/actions/workflow/status/js2264/OMA/rworkflows.yml?label=Package%20check)](https://github.com/js2264/OMA/actions/workflows/rworkflows.yml)   
📖 [Book](https://js2264.github.io/OMA/) [![deployment](https://img.shields.io/github/actions/workflow/status/js2264/OMA/pages/pages-build-deployment?label=Book%20deployment)](https://github.com/js2264/OMA/actions/workflows/pages/pages-build-deployment)  
🐳 [Docker](https://github.com/js2264/OMA/pkgs/container/OMA) [![biocbook](https://img.shields.io/github/actions/workflow/status/js2264/OMA/biocbook.yml?label=Docker%20image)](https://github.com/js2264/OMA/actions/workflows/biocbook.yml)  
<!-- badges: end -->

**README from [microbiome/OMA](https://github.com/microbiome/OMA)**

# Orchestrating Microbiome Analysis Book <img src="inst/assets/mia_logo.png" align="right" width="120" />

## Overview

This is a reference cookbook for performing **Microbiome Analysis** with 
Bioconductor in R. This is a book based on Quarto and **`BiocBook`** 
(https://www.bioconductor.org/packages/release/bioc/html/BiocBook.html).

## Deployment

The book is automatically built and deployed from the `devel` branch to 
the `gh-pages` branch using GitHub Actions.

You can also preview it locally after _cloning_ this Github
repository. This is useful for instance if you like to suggest
improvements in the material. You can use this to test the build
before making a pull request to add your new changes in the official
release.

Building and viewing the book locally involves the following steps:

1. Install the necessary dependencies to build to book, if necessary:

```
BiocManager::install(remotes::local_package_deps(dependencies=TRUE))
devtools::install('.')
```

2. Render and view the book:

```
BiocBook::preview(BiocBook::BiocBook('.'))
```

## Development and Contributions

To contribute reports, follow the Git flow procedure (you can see instructions
to [getting started with Github](https://docs.github.com/en/get-started)):

1. Fork the project
2. Clone your fork
3. Modify the material. **Check that the code runs successfully.**
4. Check locally that the changes render successfully (see above)
5. Add and commit the changes to your fork
6. Create a pull request (PR) from your fork back to the original repo
7. Fix and discuss issues in the review process

You can set OMA `devel` branch as your _upstream_ branch and pull the
changes from that before making new Pull Requests (see below). This way you can
make sure that your local version is in sync with the latest full
release.

### Setting upstream

After you forked OMA, you have two repositories to care about:

- **origin:** your own Github fork of OMA, under your github account
- **upstream:** [`devel` branch of OMA](https://github.com/microbiome/OMA/)

The origin you have after you cloned your own fork.

The upstream you can set on command line as follows, for instance (and
you can educate yourself more through various online resources on
using Git/hub):

```
git remote add upstream git@github.com:microbiome/OMA.git
```

Pull changes from the _origin_ and _upstream_ to your local version:

```
git fetch --all
git merge origin/devel
git merge upstream/devel
```

Sync your local version with the _origin_ and _upstream_:

```
git add . 
git commit -am "my changes"
```

Push your changes to origin:

```
git push origin devel
```

After this you can open a PR from origin to the [official devel branch](https://github.com/microbiome/OMA/) in Github.

### Adding new sections

Please note that chapters should be independent of each other.

- Create the relevant `.qmd` file; follow the numbering logic.
- Add it also to the list in file [inst/assets/_book.yml](inst/assets/_book.yml). 
- **Add any new dependency you use to the [DESCRIPTION](DESCRIPTION) file**.

### Styling

1. Use a coding style consistent with the rest of the book. It makes the book
look coherent.

2. Use quarto-styling instead or Rmarkdown style in code options. See
[quarto execution options](https://quarto.org/docs/computations/execution-options.html).

3. All chunks should have labels. It is easier to see for instance which chunk
is causing problems.

4. Focus on interpreting the results, as this is often what readers find most
helpful and clarifying.

5. Use callout blocks where possible. They enhance the book's appearance and
help structure the text effectively. See
[quarto callout blocks](https://quarto.org/docs/authoring/callouts.html).

6. Use the following style: \`package\` and \`function()\`.

7. The text lines should not exceed 80 characters.

8. The code lines should not exceed 80 characters. Code exceeding this limit
does not fit to html pages correctly.

Instead of this

```
long_variable_name <- long_function_name(parameter = "parameter".
                                         parameter2 = "parameter2")
```

use the following styling

```
long_variable_name <- long_function_name(
    parameter = "parameter".
    parameter2 = "parameter2"
	)
```

	

# Code of conduct

Please note that the OMA project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
