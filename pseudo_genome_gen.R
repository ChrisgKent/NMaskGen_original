#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(argparser))
suppressMessages(library(msa))


p <- arg_parser("pseudo genome Generator") 
# Adding flags
p <- add_argument(p, "--input", help="input directory")
p <- add_argument(p, "--output", help="output directory")
p <- add_argument(p, "--repair", help="Use a ref genome to repair 5' and 3' regions (TRUE/FALSE)", default = TRUE)
p <- add_argument(p, "--repair_genome", help= "The dir of the ref genome", default = "resources/MN908947.3.fasta")
p <- add_argument(p, "--m", help= "Number of Cores (Windows OS has to use 1)", default = 1)
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
if(dir.exists(argv$output)){# If the dir exists
  if(!dir.exists(paste0(argv$output, "/tmp"))){dir.create(paste0(argv$output, "/tmp"))}
}else{# if it doesn't exists
  cat(paste0("Creating output directory: ", argv$output), "\n")
  dir.create(argv$output)
  dir.create(paste0(argv$output, "/tmp"))
}

pseudo_gen <- function(x){
  #Functions 
  continuous_func <- function(x){
    group <- 1
    # Determins if its a continuous insertion, groups continuous ones
    for(i in 1:nrow(x)){
      if(i == 1){x$group[i] <- group
      }else if(x$start[i] - x$start[i-1] == 1){x$group[i] <- group
      }else{group <- group +1
      x$group[i] <- group}
    }
    x
  }
  indel_func <- function(x,y = 4){
    indel <- character()
    for(i in unique(x$group)){
      for_dat <- filter(x, group == i)
      indel[i] <- Biostrings::subseq(repair2[[1]], start = min(for_dat$start)-y, end = max(for_dat$start)+y) %>% as.character()
    }
    indel
  }
  
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
  
  rm(sub_seq_set)
  # Uses the consensus Sequence to create an unanchored pseudo Genome
  msa_consen <- msa::msaConsensusSequence(msa_seqs) %>% 
    str_replace_all("\\?", "N") %>% # If there are two bases with equal freq "?" will be returned. This is turned into an N
    str_remove_all("-") %>% # If one seq has an insert a "--" will be inserted. If the insert was in more samples the bases would be returned.
    BStringSet()
  names(msa_consen) <- paste0(x, "_pseudo_genome")
  
  if(repair){
    cat(paste0("Repairing Genome ", x, "\n"))
    # Aligns the pseudo to the referance
    psuedo_ref <- BStringSet(c(msa_consen, repair_ref))
    repair1 <- msaClustalOmega(psuedo_ref, type = "dna")
    repair2 <- BStringSet(repair1)
    unrepaired <- repair2
    
    # Determins if both sequences agree at each position
    conMat <- consensusMatrix(repair1)
    conMat <- conMat[rowSums(conMat) !=0,]
    conMax <- integer()
    for(i in 1:ncol(conMat)){
      conMax[i] <- max(conMat[,i])
    }
    repair_dat <- data.frame(pos= 1:ncol(conMat), concen = conMax)
    
    # Finds the first site of consensus between ref and pseudo
    f_2 <- repair_dat[repair_dat$concen == 2,] %>% .[,-2] %>% min()
    
    # Finds the last site of consensus between ref and pseudo
    l_2 <- repair_dat[repair_dat$concen == 2,] %>% .[,-2] %>% max()
    
    # Repairs the left end
    ## By replacing the pseudo genomes sequences with the ref upto the first site of consensus. This should just replace the tails
    repair2[[1]][1:f_2-1] <- repair2[[2]] %>% Biostrings::subseq(start = 1, end = f_2-1)
    # Repairs the end 
    repair2[[1]][(l_2+1):length(repair2[[1]])] <- repair2[[2]] %>% Biostrings::subseq(start = (l_2+1), end = length(repair2[[1]]))
    
    # Determing if there are any insertions into the pseudo, relative to the ref
    del <- matchPattern("-", repair2[[1]])
    if(length(del) != 0){# There has been an insertion
      cat(paste0("Deletion detected: ", x, "\n"))
      del_df <- as.data.frame(del@ranges) %>% mutate(group = 0)
      del_df <- continuous_func(del_df)
      
      # Saves the bases prior to repair, for use in log file
      delseq <- indel_func(del_df)
      
      # Repairs the deletions by masking all "-" with "N" with in the psudeo genome
      repair2[[1]] <- repair2[[1]] %>% str_replace_all("-", "N")
      
      # Saves the bases post repair, for use in log file
      repdelseq <- indel_func(del_df)
    }
    
    insert <- matchPattern("-", repair2[[2]])
    if(length(insert) != 0){# There has been an insertion
      cat(paste0("Insertion detected: ", x, "\n"))
      insert_df <- as.data.frame(insert@ranges) %>% mutate(group = 0)
      insert_df <- continuous_func(insert_df)
      
      # Saves the bases prior to repair, for use in log file
      buffer <- 4
      insertseq <- indel_func(insert_df,buffer)
      
      # Replaces all gaps with the * symble, then removes all to prevent any shifts in the index
      for(i in unique(insert_df$group)){
        tmp <- filter(insert_df, group == i)
        min_start = min(tmp$start)
        max_start = max(tmp$start)
        repair2[[1]][min_start:max_start] <- BString(paste0(rep("*", length(tmp$start)), collapse =""))
        repair2[[1]][min_start-1] <- BString("N")
      }
      repair2[[1]] <- str_remove_all(repair2[[1]], "\\*")
      repair2[[2]] <- str_remove_all(repair2[[2]], "-")
      
      # Saves the bases post repair, for use in log file
      repinsertseq <- character()
      for(i in 1:length(insertseq)){
        insert_str <- str_split(insertseq[i], "", simplify = TRUE)
        insert_str[buffer] <- "N"
        
        repinsertseq[i] <- insert_str[c(1:4, seq(length(insert_str)-buffer+1, length(insert_str)))] %>% 
          paste0(collapse = "")
      }
    }
    
    # Returns the repaired genome
    msa_consen_repair <- repair2[[1]] %>% BStringSet()
    names(msa_consen_repair) <- names(repair2)[1]
    
    # Writes the repaired Genome
    writeXStringSet(msa_consen_repair, paste0(argv$output, "/", x, "_pseudo_repaired.fasta"))
                                                                                                                                                                                                              
  }
  
  if(bed){
    bed_consen <- consensusMatrix(repair2)
    
    max_prop <- numeric()
    for(L in 1:ncol(bed_consen)){
      max_prop[L] <- max(bed_consen[,L])
    }
    
    data.frame(pos_1base = 1:length(repair2[[1]]),
                          dif = max_prop !=2) %>% 
      filter(dif == TRUE) %>%
      mutate(chrom = names(repair_ref),
             chromStart = pos_1base-1,
             chromEnd = pos_1base) %>% 
      select(chrom,chromStart,chromEnd) %>%
      write_delim(., 
                  delim = "\t", 
                  file = paste0(argv$output, "/", x, "_pseudo_genome.bed"),
                  col_names = FALSE)
    }
  
  # Generates and writes a log file containing all seqs that were used in pseudo Genome generation
  if(repair){
    left_msa_text <- unrepaired %>% 
      Biostrings::subseq(start =1, end = f_2+5) %>% 
      as.character() 
    
    right_msa_text <- unrepaired %>% 
      Biostrings::subseq(start = (l_2+1)-5, end = length(repair2[[1]])) %>% 
      as.character() 
    
    log <- paste(paste0("Seqs used in ", x, "_pseudo_genome:\n"), paste0(sub_seqs_names, collapse = "\n")) %>%
      paste0("\n",., "\n\nRepaired using ", names(repair_ref),
             "\n\nREPAIR: \nA 5 base overlap is shown which was not repaired \n\nLeft pseudo: \n ",
             left_msa_text[1], "\n",
             "Left ref genome: \n",
             left_msa_text[2],
             "\n\n",
             "Right pseudo: \n",
             right_msa_text[1], "\n",
             "Right ref genome: \n",
             right_msa_text[2])
    
    
    # Adding any insertions into the log file 
    if(length(insert) != 0){
      insrt_txt <- "\n\nINSERTIONS\n\nOriginal\t\t\tReplaced\t\t\tSequences"
      for(i in unique(insert_df$group)){
        for_dat <- filter(insert_df, group == i)
        
        #Find the new index locations of the inserts
        ts <- matchPattern(repinsertseq[i], repair2[[1]]) %>% 
          .@ranges %>%
          as.data.frame()
        rep_insert_log <- Biostrings::subseq(repair2, start = ts$start, end = ts$end) %>%
          as.character()
        
        #Find the index of the insert locations on the unrepaired genome
        ts <- matchPattern(insertseq[i], unrepaired[[1]]) %>% 
          .@ranges %>%
          as.data.frame()
        insert_log <- Biostrings::subseq(unrepaired, start = ts$start, end = ts$end) %>%
          as.character()
        
        insrt_txt[i+1] <- paste0("\n",insert_log[1],"\t\t\t",rep_insert_log[1],"\t\t\t",names(repair2)[1],"\n",
                                 insert_log[2],"\t\t\t",rep_insert_log[2],"\t\t\t",names(repair2)[2],"\n")
      }
     
      log <- paste0(log, paste0(insrt_txt, collapse = ""))
      
    }else{
      log <- paste0(log, "\n\n", "NO INSERTIONS")
    }
    
    # Adding any deletions into the log file 
    if(length(del) != 0){
      del_txt <- "\n\nDELETIONS\n\nOriginal\t\t\tReplaced\t\t\tSequences"
      for(i in unique(del_df$group)){
        for_dat <- filter(del_df, group == i)
        
        #Find the new index locations of the inserts
        ts <- matchPattern(repdelseq[i], repair2[[1]]) %>% 
          .@ranges %>%
          as.data.frame()
        rep_del_log <- Biostrings::subseq(repair2, start = ts$start, end = ts$end) %>%
          as.character()
        
        #Find the index of the insert locations on the unrepaired genome
        ts <- matchPattern(delseq[i], unrepaired[[1]]) %>% 
          .@ranges %>%
          as.data.frame()
        del_log <- Biostrings::subseq(unrepaired, start = ts$start, end = ts$end) %>%
          as.character()
        
        del_txt[i+1] <- paste0("\n",del_log[1],"\t\t\t",rep_del_log[1],"\t\t\t",names(repair2)[1],"\n",
                                 del_log[2],"\t\t\t",rep_del_log[2],"\t\t\t",names(repair2)[2],"\n")
      }
      
      log <- paste0(log, paste0(del_txt, collapse = ""))
      
    }else{log <- paste0(log, "\n\nNO DELETIONS")}
    
    #Final Validation of ref genome repair
    log <- paste0(log,paste0("\n\nFinal Validation of repair:\n",repair2[[2]] == repair_ref))
    
    }else{
    log <- paste(paste0("Seqs used in ", x, "_pseudo_genome: \n"), paste0(sub_seqs_names, collapse = "\n")) 
    writeXStringSet(msa_consen, paste0(argv$output, "/", x, "_pseudo.fasta"))
    }
  
  cat(paste0("Saving files for ", x),"\n")
  write_file(log,  paste0(argv$output, "/",  x, "_logfile.txt"))

}

# To ensure flexability, asking for 1 core does not use the mc.apply()
if(argv$m == 1){
  dat <- lapply(pango_list, pseudo_gen)
}else{
  dat <- mclapply(pango_list, pseudo_gen, mc.cores = as.numeric(argv$m))
}

