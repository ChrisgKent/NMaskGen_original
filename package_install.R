cran <- list("tidyverse",
          "BiocManager",
          "readxl",
          "argparser")

bioc <- list("Biostrings",
          "msa")

cran_installer <- function(x){
  if(!require(x, character.only = TRUE)){
    install.packages(x, dependencies = TRUE)
    cat(paste0("Installing: ", x, " from CRAN \n"))
  }else{cat(paste0(x, " is already installed \n"))}}



bioc_installer <- function(x){
  if(!require(x, character.only = TRUE)){
    BiocManager::install(x)
    cat(paste0("Installing: ", x, " from Bioconductor \n"))
  }else{cat(paste0(x, " is already installed \n"))}}

suppressMessages(dat <- lapply(cran, cran_installer))
suppressMessages(library(BiocManager))
suppressMessages(dat <- lapply(bioc, bioc_installer))

cat("All packages installed")

