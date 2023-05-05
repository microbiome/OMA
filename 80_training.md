# (PART) Training {-}

# Training {#training}

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

The page provides practical information to support training and self-study.


## Checklist {#checklist}

Brief checklist to prepare for training (see below for links).

 * Install the recommended software
 * Watch the short online videos and familiarize with the other available material
 * Join Gitter online chat for support


## Recommended software {#software}

We recommend to install and set up the relevant software packages on
your own computer as this will support later use. The essential
components to install include:

* [R (the latest official release)](https://www.r-project.org/) 

* [RStudio](https://www.rstudio.com/products/rstudio/download/);
  choose "Rstudio Desktop" to download the latest version. Check the
  [Rstudio home page](https://www.rstudio.com/) for more
  information. RStudio is optional.

* Install key R packages (Section \@ref(packages) provides an installation script)


* After a successful installation you can consider trying out examples
  from Section \@ref(exercises) already before training. **You can run
  the workflows by simply copy-pasting examples.** You can then test
  further examples from this tutorial, modifying and applying these
  techniques to your own data. Plain source code for the individual chapters of this book are available via [Github](https://github.com/microbiome/OMA/tree/master/R)

## Study material {#material}

We encourage to familiarize with the material and test examples in advance.

 * [Short online videos](https://www.youtube.com/playlist?list=PLjiXAZO27elAJEptP59BN3whVJ61XIkST) on microbiome data science with R/Bioconductor  

 * [Lecture slides](https://microbiome.github.io/)

 * [Orchestrating Microbiome Analysis with R/Bioconductor (OMA)](https://microbiome.github.io/OMA/) (this book)

 * [Exercises](https://microbiome.github.io/OMA/exercises.html) for self-study

 * [Resources](https://microbiome.github.io/OMA/resources.html) and links to complementary external material



## Support and resources

For online support on installation and other matters, join us at
[Gitter](https://gitter.im/microbiome/miaverse?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge).

You are also welcome to connect through various channels with our
broader [developer and user community](https://microbiome.github.io).


## Code of Conduct {#coc}

We support the [Bioconductor Code of Conduct](https://bioconductor.github.io/bioc_coc_multilingual/). The community values an open approach to science that promotes 

  - sharing of ideas, code, software and expertise
  - a kind and welcoming environment, diversity and inclusivity
  - community contributions and collaboration
