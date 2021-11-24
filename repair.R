repair_func <- function(name, seq, ref, log = FALSE, log_dir, sub_seqs_names, base_check){
  #
  # This function takes a genome (seq) and a referance (ref)
  # It returns an altered seq, so it has the same index system as the ref.
  # 
  # Seq and ref need to be BStringSets (or DNA) with length 1. So use my_set[1] rather than my_set[[1]]
  #
  #
  # This is done by;
  #   Duplicating the 5' and 3' ends of the ref seq onto the seq, only were there is no consensus.
  #         seq:  ----AGT
  #         ref:  AGCTAGT
  #       repair: AGCTAGT
  #   Masking deletions locations with Ns
  #         seq:  ACGT------GTCA
  #         ref:  ACGTGTACGAGTCA
  #       repair: ACGTNNNNNNGTCA
  #   Masking insertions by removing the insert and masking the -1 base with an N
  #         seq:  ACGTGTACGAGTCA
  #         ref:  ACGT------GTCA
  #       repair: ACGNGTCA
  #
  

  cat(paste0("Repairing Genome ", name, "\n"))
    
  # Aligns the pseudo to the referance
    psuedo_ref <- BStringSet(c(seq, ref))
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
    
    if(consensusMatrix(subseq(repair2, start = 1, width = 1))[1,1] != 2){
      # Finds the first site of consensus between ref and pseudo
      f_2 <- repair_dat[repair_dat$concen == 2,] %>% .[,-2] %>% min()
      # Repairs the right end 
      repair2[[1]][1:f_2-1] <- repair2[[2]] %>% Biostrings::subseq(start = 1, end = f_2-1)
    }
    
    if(consensusMatrix(subseq(repair2,end = length(repair2[[2]]), width = 1))[1,1] != 2){
      # Finds the last site of consensus between ref and pseudo
      l_2 <- repair_dat[repair_dat$concen == 2,] %>% .[,-2] %>% max()
      # Repairs the left end 
      repair2[[1]][(l_2+1):length(repair2[[1]])] <- repair2[[2]] %>% Biostrings::subseq(start = (l_2+1), end = length(repair2[[1]]))
    }
    
    # Determing if there are any insertions into the pseudo, relative to the ref
    del <- matchPattern("-", repair2[[1]])
    if(length(del) != 0){# There has been an insertion
      cat(paste0("Deletion detected: ", name, "\n"))
      del_df <- as.data.frame(del@ranges) %>% mutate(group = 0)
      del_df <- continuous_func(del_df)
      
      # Saves the bases prior to repair, for use in log file
      delseq <- indel_func(del_df, seq = repair2[[1]])
      
      # Repairs the deletions by masking all "-" with "N" with in the psudeo genome
      repair2[[1]] <- repair2[[1]] %>% str_replace_all("-", "N")
      
      # Saves the bases post repair, for use in log file
      repdelseq <- indel_func(del_df, seq = repair2[[1]])
    }
    
    insert <- matchPattern("-", repair2[[2]])
    if(length(insert) != 0){# There has been an insertion
      cat(paste0("Insertion detected: ", name, "\n"))
      insert_df <- as.data.frame(insert@ranges) %>% mutate(group = 0)
      insert_df <- continuous_func(insert_df)
      
      # Saves the bases prior to repair, for use in log file
      buffer <- 4
      insertseq <- indel_func(insert_df,buffer, seq = repair2[[1]])
      
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
    msa_consen_repair <- repair2 %>% BStringSet()
    names(msa_consen_repair) <- names(repair2)
    
    if(log){
      left_msa_text <- unrepaired %>% 
        Biostrings::subseq(start =1, end = f_2+5) %>% 
        as.character() 
      
      right_msa_text <- unrepaired %>% 
        Biostrings::subseq(start = (l_2+1)-5, end = length(repair2[[1]])) %>% 
        as.character() 
      
      log <- paste(paste0("Seqs used in ", name, "_pseudoref:\n"), paste0(sub_seqs_names, collapse = "\n")) %>%
        paste0("\n",., "\n\nRepaired using ", names(ref),
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
      log <- paste0(log,paste0("\n\nFinal Validation of repair:\n",repair2[[2]] == ref))
      
      if(sum(base_check[base_check$names == "consen",-1]) != 0){
        log <- paste0(log, "\n\nNon [ACGT] bases found in ", name, " consensus sequence\n\n")
      }
      
      write_file(log,  paste0(log_dir, "/",  name, "_logfile.txt"))
      
    }
    
    return(msa_consen_repair)
}


# These two functions are required for repair to work
continuous_func <- function(x){
  # Determins if its a continuous insertion, and groups continuous ones
  group <- 1
  for(i in 1:nrow(x)){
    if(i == 1){x$group[i] <- group
    }else if(x$start[i] - x$start[i-1] == 1){x$group[i] <- group
    }else{group <- group +1
    x$group[i] <- group}
  }
  return(x)
}

indel_func <- function(x,y = 4, seq = repair2[[1]]){
  # Saves the nuclotide seq of the indel for later
  indel <- character()
  for(i in unique(x$group)){
    for_dat <- filter(x, group == i)
    indel[i] <- Biostrings::subseq(seq, start = min(for_dat$start)-y, end = max(for_dat$start)+y) %>% as.character()
  }
  indel
}