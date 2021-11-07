#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(argparser))
suppressMessages(library(msa))


p <- arg_parser("pseudo genome Generator") 
# Adding flags
p <- add_argument(p, "--input", help="input directory")
p <- add_argument(p, "--output", help="output directory")
p <- add_argument(p, "--repair", help="Use a ref genome to repair 5' and 3' regions (TRUE/FALSE)")
p <- add_argument(p, "--repair_genome", help= "The dir of the ref genome")
p <- add_argument(p, "--mc", help= "Number of Cores")

argv <- parse_args(p)

# If the repair argument is TRUE reads in the referance genome wanted 
repair <- as.logical(argv$repair)
if(repair){
  repair_ref <- readBStringSet(argv$repair_genome) %>%
    DNAStringSet()
  cat("Reading in Ref Genome")
  }


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
  
  if(repair){ 
    cat(paste0("Repairing Genome ", x, "\n"))
    # Aligns the pseudo to the referance
    psuedo_ref <- BStringSet(c(msa_consen, repair_ref))
    repair1 <- msaClustalOmega(psuedo_ref, type = "dna")
    repair2 <- BStringSet(repair1)
    
    # Determins if both seqeucnes agree at each position
    conMat <- consensusMatrix(repair1)
    conMat <- conMat[rowSums(conMat) !=0,]
    conMax <- integer()
    for(i in 1:ncol(conMat)){
      conMax[i] <- max(conMat[,i])
    }
    repair_dat <- data.frame(pos= 1:ncol(conMat), concen = conMax)
    
    # Finds the first site of consensus between ref and pysudo
    f_2 <- repair_dat[repair_dat$concen == 2,] %>% .[,-2] %>% min()
    
    # Finds the last site of Consen between ref and pysudo
    l_2 <- repair_dat[repair_dat$concen == 2,] %>% .[,-2] %>% max()
    
    # Repairs the start
    ## By replacing the pseudo genomes sequecnes with the ref upto the first site of consensus. This should just replace the tails
    repair2[[1]][1:f_2-1] <- repair2[[2]] %>% Biostrings::subseq(start = 1, end = f_2-1)
    # Repairs the end 
    repair2[[1]][(l_2+1):length(repair2[[1]])] <- repair2[[2]] %>% Biostrings::subseq(start = (l_2+1), end = length(repair2[[1]]))
    cat(paste0("Returning Repaired ", x, "\n"))
    # Returns the repaired genome
    msa_consen_repair <- repair2[[1]] %>% str_remove_all("-") %>% BStringSet()
    names(msa_consen_repair) <- names(repair2)[1]
    
    # Writes the repaired Genome
    writeXStringSet(msa_consen_repair, paste0(argv$output, "/", x, "_pseudo_repaired.fasta"))
                                                                                                                                                                                                              
  }
  

  # Generates and writes a log file containing all seqs that were used in pseudo Genome generation
  if(repair){
    log <- paste(paste0("Seqs used in ", x, "_pseudo_genome: \n"), paste0(sub_seqs_names, collapse = "\n")) %>%
      paste0(., "\n Repaired using ", names(repair_ref))
    }else{
    log <- paste(paste0("Seqs used in ", x, "_pseudo_genome: \n"), paste0(sub_seqs_names, collapse = "\n")) 
    }
  cat(paste0("Saving files for ", x),"\n")
  writeXStringSet(msa_consen, paste0(argv$output, "/", x, "_pseudo.fasta"))

  write_file(log,  paste0(argv$output, "/",  x, "_logfile.txt"))

}


dat <- mclapply(pango_list, pseudo_gen, mc.cores = as.numeric(argv$mc))
