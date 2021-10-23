#!/usr/bin/env Rscript
library(tidyverse)
library(Biostrings)

files <- list.files("UK_scheme_data/") %>% 
  .[str_detect(., ".csv") == FALSE] %>%
  .[str_detect(., ".xlsx") == FALSE]

pango_list <- split(files, seq(length(files)))

psesdo_gen <- function(x){
  sub_seqs_files <- list.files(paste0("UK_scheme_data/", x))
  sub_seqs_names <- sub_seqs_files %>% str_remove(".fasta")
  
  sub_seq_set <- DNAStringSet()
  for(i in 1:length(sub_seqs_files)){
    sub_seq_set[i] <- readBStringSet(paste0("UK_scheme_data/", x, "/", sub_seqs_files[i]))
  }
  names(sub_seq_set) <- sub_seqs_names
  
  msa_seqs <- msa::msaClustalOmega(sub_seq_set)
  msa_consen <- msa::msaConsensusSequence(msa_seqs) %>% DNAStringSet()
  names(msa_consen) <- paste0(x, "_pseudo_genome")
  
  writeXStringSet(msa_consen, paste0("UK_scheme_data/pseudo_genomes/", x, "_pseudo_genome.fasta"))
  log <- paste(paste0("Seqs used in ", x, "_pseudo_genome: \n"), paste0(sub_seqs_names, collapse = "\n"))
  
  write_file(log,  paste0("UK_scheme_data/pseudo_genomes/", x, "_logfile.txt"))

}

mclapply(pango_list, psesdo_gen)




