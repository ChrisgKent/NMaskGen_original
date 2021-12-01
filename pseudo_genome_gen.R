#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(argparser))
suppressMessages(library(msa))


p <- arg_parser("Pseudo Genome Generator",hide.opts = TRUE) 
# Adding flags
p <- add_argument(p, "--input", help="input directory")
p <- add_argument(p, "--output", help="output directory")
p <- add_argument(p, "--repair", help="Use a ref genome to repair 5' and 3' regions (TRUE/FALSE)", default = TRUE)
p <- add_argument(p, "--repair_genome", help= "The dir of the ref genome", default = "resources/MN908947.3.fasta")
p <- add_argument(p, "--cores", help= "Number of Cores (Windows OS has to use 1)", default = 1)
p <- add_argument(p, "--bed", help= "Should a .bed file written (requires --repair = TRUE)", default = TRUE)

argv <- parse_args(p)

bed <- as.logical(argv$bed)
# If the repair argument is TRUE reads in the referance genome wanted 
repair <- as.logical(argv$repair)
if(repair){
  repair_ref <- readBStringSet(argv$repair_genome) %>%
    DNAStringSet()
  cat("Reading in Ref Genome \n")
  }

# Finds all PANGO lineages dirs
files <- list.files(argv$input) %>% 
  .[str_detect(., ".csv") == FALSE] %>%
  .[str_detect(., ".xlsx") == FALSE] %>% 
  .[str_detect(., "pseudo_genomes") == FALSE]
  
## Generates them into a list for lapply()
pango_list <- split(files, seq(length(files)))

# Checks if the output dir already exists, if not it generates it and a tmp directory 
if(!dir.exists(argv$output)){# If the dir doesn't exist
  cat(paste0("Creating output directory: ", argv$output), "\n")
  dir.create(argv$output)
  dir.create(paste0(argv$output, "/tmp"))
}else{
  if(!dir.exists(paste0(argv$output, "/tmp"))){
    dir.create(paste0(argv$output, "/tmp"))
  }
}

pseudo_gen <- function(x){
  # Finds all seqs for each PANGO and ensures there's more than one.
  sub_seqs_files <- list.files(paste0(argv$input,"/", x))
  if(length(sub_seqs_files) <= 1){stop(paste0("Not enough Sequences detected", x))}
  sub_seqs_names <- sub_seqs_files %>% str_remove(".fasta")
  
  # Reads in the sequences as a DNA string object
  sub_seq_set <- DNAStringSet()
  for(i in 1:length(sub_seqs_files)){
    sub_seq_set[i] <- readBStringSet(paste0(argv$input,"/", x, "/", sub_seqs_files[i]))
  }
  names(sub_seq_set) <- sub_seqs_names
  
  # Runs MSA clustal Omega on the sequences
  cat(paste0("Running ClustalOmega on ", x), "\n")
  msa_seqs <- msa::msaClustalOmega(sub_seq_set, type = "dna")
  
  writeXStringSet(BStringSet(msa_seqs), paste0(argv$output, "/tmp/", x, "_msa.fasta"))
  
  # Uses the consensus Sequence to create an unanchored pseudo Genome
  source("consensus_gen.R")
  msa_consen <- consen_gen(DNAStringSet(msa_seqs)) 
  names(msa_consen) <- paste0(x, "_pseudoref")
  
  # Checks for non standard bases
  base_check <- Biostrings::letterFrequency(sub_seq_set, DNA_ALPHABET) %>% .[,-c(1:4)]
  base_check <- rbind(base_check, Biostrings::letterFrequency(msa_consen, DNA_ALPHABET) %>% .[,-c(1:4)]) %>% as.data.frame()
  base_check <- cbind(names = c(sub_seqs_names, "consen"), base_check) %>% as.data.frame()
  write_csv(base_check, paste0(argv$output, "/tmp/", x, "_basecheck.csv"))
  
  if(sum(base_check[base_check$names == "consen",-1]) != 0){
    cat(paste0("Non [ACGT] bases found in ", x, " consensus sequence"))
  }
  
  if(repair){
    source("repair.R")
    repair2 <- repair_func(seq = msa_consen, 
                           ref = repair_ref, 
                           name = x, 
                           log = TRUE, 
                           log_dir = argv$output, 
                           sub_seqs_names = sub_seqs_names,
                           base_check = base_check)
    
    writeXStringSet(repair2[1], paste0(argv$output, "/", x, "_pseudoref.fasta"))
  }else{
    log <- paste(paste0("Seqs used in ", x, "_pseudoref: \n"), paste0(sub_seqs_names, collapse = "\n")) 
    writeXStringSet(msa_consen, paste0(argv$output, "/", x, "_pseudo.fasta"))
    write_file(log,  paste0(argv$output, "/",  x, "_logfile.txt"))
  }
  
  if(bed){# This section writes a bed file
    ## that contains the locations and corresponding change between the ref and the pseudo genome 
    ## Therefore, it can be used to de novo generate the pseudo genomefrom the ref and bed files
    bed_consen <- consensusMatrix(repair2)
    max_prop <- numeric()
    for(L in 1:ncol(bed_consen)){
      # Finds the most common base at each position
      max_prop[L] <- max(bed_consen[,L])
    } 
    
    # Filters for locations with variation between the two genomes 
    bed_data <- data.frame(pos_1base = 1:length(repair2[[1]]),dif = max_prop !=2) %>% 
      filter(dif == TRUE) %>%
      mutate(chrom = names(repair_ref),
             chromStart = pos_1base-1,
             chromEnd = pos_1base) %>% 
      select(chrom,chromStart,chromEnd)
    
    # For each location with variation, the ref base and pseudo base are found 
    bed_names <- character()
    for(i in 1:nrow(bed_data)){
      for_dat <- Biostrings::subseq(repair2, start = bed_data$chromStart[i]+1, width = 1)
      ref <- for_dat[[2]] %>% as.character()
      psuedo <- for_dat[[1]] %>% as.character()
      
      # The two bases are parsed into a string
      bed_names[i] <- paste0(ref, "to", psuedo)
    }
    
    # bedfile is written
    bed_data <- mutate(bed_data,name = bed_names)
    write_delim(bed_data, 
                delim = "\t", 
                file = paste0(argv$output, "/", x, "_pseudoref.bed"),
                col_names = FALSE)
    }
  
  cat(paste0("Saving files for ", x),"\n")


}

# To ensure flexability, asking for 1 core does not use the mc.apply()
if(argv$m == 1){
  dat <- lapply(pango_list, pseudo_gen)
}else{
  dat <- parallel::mclapply(pango_list, pseudo_gen, mc.cores = as.numeric(argv$m))
}


