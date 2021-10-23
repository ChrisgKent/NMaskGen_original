#!/usr/bin/env Rscript
library(tidyverse)
library(Biostrings)
library(argparser)

p <- arg_parser("pseudo genome Generator") 
# Adding flags
p <- add_argument(p, "--input", help="input directory")
p <- add_argument(p, "--output", help="output directory")

argv <- parse_args(p)

# Finds all PANGO lineages dirs
files <- list.files(argv$input) %>% 
  .[str_detect(., ".csv") == FALSE] %>%
  .[str_detect(., ".xlsx") == FALSE] %>% 
  .[str_detect(., "pseudo_genomes") == FALSE]
  
## Generates them into a list for lapply()
pango_list <- split(files, seq(length(files)))


pseudo_gen <- function(x){
  # Finds all seqs for each PANGO and ensures there's more than one.
  sub_seqs_files <- list.files(paste0(argv$input,"/", x))
  if(length(sub_seqs_files) <= 1){stop(paste0("Not enough Seqs", x))}
  sub_seqs_names <- sub_seqs_files %>% str_remove(".fasta")
  
  # Stores these as a DNA string object
  sub_seq_set <- DNAStringSet()
  for(i in 1:length(sub_seqs_files)){
    sub_seq_set[i] <- readBStringSet(paste0(argv$input,"/", x, "/", sub_seqs_files[i]))
  }
  names(sub_seq_set) <- sub_seqs_names
  
  # Runs MSA clustal Omega on the sequences
  cat(paste0("Running ClustalOmega on ", x), "\n")
  msa_seqs <- msa::msaClustalOmega(sub_seq_set, type = "dna")
  
  # Uses the consensus Sequence as the pseudo Genome
  msa_consen <- msa::msaConsensusSequence(msa_seqs) %>% 
    str_replace_all("\\?", "N") %>%
    str_replace_all("-", "N") %>%
    BStringSet()
  names(msa_consen) <- paste0(x, "_pseudo_genome")
  
  # Writes the pseudo Genome as a fasta
  cat(paste0("Saving files for ", x),"\n")
  writeXStringSet(msa_consen, paste0(argv$output, "/", x, "_pseudo.fasta"))
  
  # Generates and writes a log file containing all seqs that were used in pseudo Genome generation
  log <- paste(paste0("Seqs used in ", x, "_pseudo_genome: \n"), paste0(sub_seqs_names, collapse = "\n"))
  write_file(log,  paste0(argv$output, "/",  x, "_logfile.txt"))

}


dat <- mclapply(pango_list, pseudo_gen)




