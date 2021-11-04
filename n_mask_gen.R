#!/usr/bin/env Rscript

# CRAN packages 
suppressMessages(library(argparser))
suppressMessages(library(tidyverse))

# Bioconductor Packages
suppressMessages(library(msa))
suppressMessages(library(Biostrings))

# Adding a positional argument
p <- arg_parser("N Mask Generator") 
# Adding flags
p <- add_argument(p, "--seq", help="Directory of input file containing fasta seq/s")
p <- add_argument(p, "--threshold", help="The cut off proportion of the most common base in each position")
p <- add_argument(p, "--name_output", help="The prefix of the output files")
p <- add_argument(p, "--output_dir", help="The directory of the output files")
# Adding an option to include a differant N-mask base
p <- add_argument(p, "--mask_base", help="The directory for base of the N-mask. Default is ClusterOmega Consensus")

# Adds the ability to mimic the VCF pipeline
p <- add_argument(p, "--vcf_mimic", help="If bases prior to poly-N tract should be masked. Leave blank for no")

argv <- parse_args(p)


seq_cat <- argv$seq
threshold <- argv$threshold
output_name <- argv$name_output
output_dir <- argv$output_dir

vcf <- as.logical(argv$vcf_mimic)

# Generating the output Dirs
bed_file_name <- paste0(output_dir, "/", output_name, ".bed")
consen_file_name <- paste0(output_dir, "/", output_name, "_base_mask.fasta")
n_mask_file_name <- paste0(output_dir, "/", output_name, "_n_mask.fasta")

# Reading in the --seq file
data <- readBStringSet(seq_cat, format = "fasta") 
cat(paste0(".fasta has been found and read \n"))

# Runs ClustalOmega on the input seqs
cat(paste0("Running ClustalOmega \n"))
msa_data <- msaClustalOmega(data, type = "dna")

# Generating a consensus matrix, and removes any rows that are empty
conMat <- consensusMatrix(msa_data)
conMat <- conMat[rowSums(conMat) !=0,]

# Generates the proportion of the most common base at each position
number_of_seqs <- length(msa_data@unmasked)
conMat_prop <- conMat / number_of_seqs

max_prop <- integer()
for(i in 1:ncol(conMat_prop)){
  max_prop[i] <- max(conMat_prop[,i])
}

# Generates a .bed file which contains the all locations.
bed_data <- data.frame(pos_1base = (1:length(max_prop)),
                       pos_0base = (1:length(max_prop))-1,
                       max_prop,
                       mask =  max_prop < threshold,
                       seq_source = "NA")

# Mimicing the VCF output
if(vcf){
  cat("Mimicing .VCF")
  for(i in 1:nrow(bed_data)){
    if(i > 1 & bed_data$mask[i] == TRUE & bed_data$mask[i+1] == TRUE){
      bed_data$mask[i-1] <- TRUE
      bed_data$seq_source[i-1] <- "VCF"
      }}}else{cat("Not Mimicing VCF")}

# Determins which sequence each N comes from
bed_test <- bed_data %>% 
  filter(mask == TRUE) %>%
  select(-pos_0base) 

## Turns the MSA object into a BString Set for ease of use
msa_test <- msa_data %>% BStringSet()

## For each position to be masked, the most common base is determined. 
## Any sequence that does not the most common base at the position is collapsed into a single string and added to the seq_source
suppressWarnings(
  for(i in 1:nrow(bed_test)){
    for_loop_dat <- Biostrings::subseq(msa_test, start =  bed_test$pos_1base[i], width = 1)
    for_loop_consen <- consensusMatrix(for_loop_dat)
    most_pop_symb <- rownames(for_loop_consen)[for_loop_consen == max(for_loop_consen)]
    
    concen <- paste0("consen(",most_pop_symb,")")
    
    if(bed_test$max_prop[i] == 1){
      bed_test$seq_source[i] <- paste0(concen, " | VCF")
    }else{
      res <- data.frame(names = names(for_loop_dat)[for_loop_dat != most_pop_symb],
                      base = "NA")
      for(j in 1:nrow(res)){
        res$base[j] <- for_loop_dat[names(for_loop_dat) == res$names[j]]
        }
      text <- res %>% 
        mutate(txt = paste0(names,"(", base,")")) %>%
        .$txt %>%
        unlist(use.names = FALSE) %>%
        paste0(., collapse = " | ") %>%
        paste0(concen, " | ", .)
      
      bed_test$seq_source[i] <- text}
    
})




# Determines if user has provided custom mask base.
if(is.na(argv$mask_base)){
  # If the mask_base if left empty
  consen_seq <- msa::msaConsensusSequence(msa_data) %>%
    str_replace_all("\\?", "N") %>%
    str_replace_all("-", "N") %>%
    DNAStringSet()
  names(consen_seq) <- c("ClustalOmega_consensus")
  cat(paste0("Saving read in N-Mask base \n"))
}else{
  # If the mask_base has a directory
  consen_seq <- readDNAStringSet(argv$mask_base)
  cat(paste0("Saving consensus N-Mask base \n"))
}

cat(paste0("Saving N-Mask base \n"))

# Writes the mask base file
writeXStringSet(consen_seq, consen_file_name, format = "fasta")

# Generates anf writes the bed file, using only the positions that need masking
bed_file <- bed_test %>% 
  mutate(chrom = names(consen_seq),
         chromStart = pos_1base-1,
         chromEnd = pos_1base,
         name = seq_source) %>%
  select(chrom,chromStart,chromEnd,name)

cat(paste0("Saving .bed file \n"))

write_delim(bed_file, 
            delim = "\t", 
            file = bed_file_name,
            col_names = FALSE)


# Calls bedtools to generate the N masked file
system(paste0("bedtools maskfasta -fi ", consen_file_name, " -bed ", bed_file_name, " -fo ",  n_mask_file_name))
cat(paste0("Saving N-masked .fasta \n"))



