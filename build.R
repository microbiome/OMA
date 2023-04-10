# Empty the R folder and generate new R files from Rmd
system("rm R/*")
fs <- list.files(pattern=".Rmd"); 
rm.fs <- c("tmp.Rmd", "misc.Rmd")
fs <- setdiff(fs, rm.fs)
for (f in fs) {knitr::purl(f, output = gsub("Rmd$", "R", paste0("R/", f)))}

# R files
fsr <- list.files(pattern=".R$", path="R", full.names=TRUE);
#fsr <- setdiff(fsr, c("R/dada2_workflow.R", "R/40_machine_learning.R")) # Fix later
rm.fs <- c("R/dada2_workflow.R", "R/40_machine_learning.R")
fsr <- setdiff(fsr, rm.fs) # Fix later
print("-----------")
print("Not converting the following packages:")
print(rm.fs)
print("-----------")
# Run the R files and load all libraries
# .. run them in random order to add testing power
for (f in sample(fsr)) {print(f); suppressMessages(source(f))}

# List the session libraries
# prioritize the mia packages to ensure their priority in NAMESPACE
pkgs <- (.packages())
priority <- c("mia", "miaViz", "miaTime", "miaSim")
pkgs <- c(priority, setdiff(pkgs, priority))

# Store the library list
write.table(pkgs, file = "oma_packages.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Build the book
bookdown::render_book("index.Rmd", "bookdown::gitbook")

