# The code adapted from Huber et al. 2023 
# https://www.huber.embl.de/msmb/install_packages.R


# Source location for the up-to-date package list:
packages <- url("https://raw.githubusercontent.com/microbiome/OMA/devel/oma_packages/oma_packages.csv")
#packages <- "oma_packages.csv"
# -------------------------------------------------------------------------------------------

options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
Sys.setenv(R_REMOTES_UPGRADE = "never")

pkg_type <- switch (Sys.info()["sysname"],
                    "Linux" = "source", 
                    "both")

## Function to install packages one at a time with indication of time left
## Overall probably slower than install.packages if everything works
## but doesn't require downloading all packages first before trying to install any
installer_with_progress <- function(pkgs) {
    
    if(length(pkgs) == 0) { invisible(return(NULL)) }
    
    toInstall <- pkgs
    bp <- progress::progress_bar$new(total = length(toInstall),
                                     format = "Installed :current of :total (:percent ) - current package: :package",
                                     show_after = 0,
                                     clear = FALSE)
    
    length_prev <- length(toInstall)
    fail <- NULL
    while(length(toInstall)) {
        pkg <- toInstall[1]
        bp$tick(length_prev - length(toInstall),  tokens = list(package = pkg))
        length_prev <- length(toInstall)
        tryCatch(
            suppressMessages( BiocManager::install(pkg, quiet = TRUE, update = FALSE, ask = FALSE, type = "binary") ),
            error = function(e) { fail <<- c(fail, pkg) },
            warning = function(w) { fail <<- c(fail, pkg) },
            ## remove current package, otherwise we loop in event of failure
            ## update the list to reflect any dependencies that are now installed
            finally = { toInstall <- setdiff(toInstall, installed.packages()[, "Package"]) }
        )
    }
    bp$tick(length_prev - length(toInstall),  tokens = list(package = "DONE!"))
    
    return(fail)
}

## these packages are needed prior to the installation
if(!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(c('BiocManager'), repos = "https://cloud.r-project.org",
                     quiet = TRUE, update = FALSE, ask = FALSE, type = pkg_type)
}
## update any existing packages
BiocManager::install(update = TRUE, ask = FALSE)

if(!requireNamespace("remotes", quietly = TRUE)) {
    install.packages(c('remotes'), quiet = TRUE, update = FALSE, ask = FALSE, type = pkg_type)
}
if(!requireNamespace("magrittr", quietly = TRUE)) {
    BiocManager::install('magrittr', quiet = TRUE, update = FALSE, ask = FALSE, type = pkg_type)
}

# ---------------------------

## list of packages required for each chapters
pkgs_all <- read.table(packages)[,1]

# This will be installed manually later in this script
# since it requires the latest devel update
pkgs_all <- setdiff(pkgs_all, "Maaslin2") 

# Customization
# Github packages must be installed separately
pkgs_github <- c("miaTime", "ggord")
pkgs_nongithub <- setdiff(pkgs_all, pkgs_github)

# Maaslin2 needs an update, see
# https://forum.biobakery.org/t/xtfrm-error-with-maaslin2-default-example-in-r/5216/3
remotes::install_github("biobakery/Maaslin2")

# ---------------------------

# pkgs <- readRDS("oma_packages.rds") # Just do all at once
chapter_pkgs <- list(all=pkgs_all)
# Can check later how packages can be splitted by chapter
#chapter_pkgs <- split(pkgs$packages, pkgs$chapter)
### subset a selection of chapters if specified
#if(exists('chapter_index') && is.numeric(chapter_index)) {
#  chapter_pkgs <- chapter_pkgs[ chapter_index ]
#}

for(i in seq_along(chapter_pkgs)) {
    message("### CHAPTER: ", i, " ###")
    pkgsAvailable <- installed.packages()[, "Package"]
    pkgsToInstall <- setdiff(chapter_pkgs[[i]], c(pkgsAvailable, pkgs_github))
    BiocManager::install(pkgsToInstall, update = FALSE, upgrade = FALSE, ask = FALSE, type = pkg_type)
}

# Github packages
devtools::install_github("microbiome/miaTime")
devtools::install_github("fawda123/ggord")

## report packages no installed
## find only those not currently installed
pkgsAvailable <- installed.packages()[, "Package"]
pkgsNeeded <- unique(unlist(chapter_pkgs))
pkgsToInstall <- setdiff(pkgsNeeded, pkgsAvailable)
if(length(pkgsToInstall)) {
    message("The following packages failed to install: \n",
            paste(pkgsToInstall, collapse = ", "))
    message("You can try re-running this installation script.\n",
            "It will only try to install the missing packages.\n",
            "This may make it easier to see the information R gives about why the installation failed.\n",
            "Please contact microbiome@utu.fi for additional help.")
}

Sys.unsetenv("R_REMOTES_UPGRADE")
