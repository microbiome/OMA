
# the two libraries you need :
# library(quarto)
# library(knitr)

#all the files in the _bookdown.yml without subject_geography.Rmd
rmd_files <- c(
  "index.Rmd", "01_intro.Rmd.Rmd", "04_containers.Rmd",
  "06_packages.Rmd", "OMA.Rmd", "10_manipulation.Rmd",
  "11_taxonomic_information.Rmd", "12_quality_control.Rmd", "14_alpha_diversity.Rmd",
  "19_visualization_techniques.Rmd", "publication_interval.Rmd", "20_beta_diversity.Rmd",
  "21_microbiome_community.Rmd", "22_communitytyping.Rmd", "assay_analyses.Rmd",
  "24_biclustering.Rmd", "30_differential_abundance.Rmd", "40_machine_learning.Rmd",
  "80_training.Rmd","90_acknowledgments.Rmd","95_resources.Rmd","96_function_reference.Rmd",
  "97_extra_materials.Rmd","98_exercises.Rmd","99_bibliography.Rmd","add-comm-typing.Rmd",
  "24_clustering.Rmd","10_manipulation.Rmd","dada2_workflow.Rmd","todo.Rmd","Session_info.Rmd"
)



#The following loop wil convert the Rmd files into quarto
for (file in rmd_files) {
  # name for the qmd files : we just change the file's extension
  qmd_file <- gsub("\\.Rmd$", ".qmd", file) 

  tryCatch({
    # convert the .Rmd file into quarto format, the output is just the name of the new .qmd file created
    knitr::convert_chunk_header(file, output = qmd_file)
    cat("created:", qmd_file, "\n")
  }, error = function(e) {
    cat("Failed to convert:", file, "\n")
  })
}

# all the Rmd files int the current dir
qmd_files <- list.files(pattern = "\\.qmd$")

# de-comment the following line to render the document in the current directory
quarto::quarto_render() 