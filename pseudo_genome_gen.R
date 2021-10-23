#!/usr/bin/env Rscript
library(tidyverse)
library(Biostrings)

# Finds all PANGO lineages dirs
files <- list.files("UK_scheme_data/") %>% 
  .[str_detect(., ".csv") == FALSE] %>%
  .[str_detect(., ".xlsx") == FALSE] %>% 
  .[str_detect(., "pseudo_genomes") == FALSE]
  
## Generates them into a list for lapply()
pango_list <- split(files, seq(length(files)))


pseudo_gen <- function(x){
  # Finds all seqs for each PANGO and ensures there's more than one.
  sub_seqs_files <- list.files(paste0("UK_scheme_data/", x))
  if(length(sub_seqs_files) <= 1){stop(paste0("Not enough Seqs", x))}
  sub_seqs_names <- sub_seqs_files %>% str_remove(".fasta")
  
  # Stores these as a DNA string object
  sub_seq_set <- DNAStringSet()
  for(i in 1:length(sub_seqs_files)){
    sub_seq_set[i] <- readBStringSet(paste0("UK_scheme_data/", x, "/", sub_seqs_files[i]))
  }
  names(sub_seq_set) <- sub_seqs_names
  
  # Runs MSA clustal Omega on the sequences
  cat(paste0("Running ClustalOmega on ", x))
  msa_seqs <- msa::msaClustalOmega(sub_seq_set)
  
  # Uses the consensus Sequence as the pseudo Genome
  msa_consen <- msa::msaConsensusSequence(msa_seqs) %>% BStringSet()
  names(msa_consen) <- paste0(x, "_pseudo_genome")
  
  # Writes the pseudo Genome as a fasta
  cat(paste0("Saving files for ", x))
  writeXStringSet(msa_consen, paste0("UK_scheme_data/pseudo_genomes/", x, "_pseudo_genome.fasta"))
  
  # Generates and writes a log file containing all seqs that were used in pseudo Genome generation
  log <- paste(paste0("Seqs used in ", x, "_pseudo_genome: \n"), paste0(sub_seqs_names, collapse = "\n"))
  write_file(log,  paste0("UK_scheme_data/pseudo_genomes/", x, "_logfile.txt"))

}

mclapply(pango_list, pseudo_gen)




