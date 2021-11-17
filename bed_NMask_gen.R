#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(argparser))

# This can generate the psuedo genome from the bed file produced via pseudo_genome_gen.R
## Currently being used as a sanity checker

p <- arg_parser("De novo pseudo genome generator") 
# Adding flags
p <- add_argument(p, "--ref", help="Directory of the reference genome",  default = "resources/MN908947.3.fasta")
p <- add_argument(p, "--bed", help= "Directory that contains the .bed files")
p <- add_argument(p, "--output", help= "Directory that output files are written", default = "NMask")
p <- add_argument(p, "--name", help= "Output file prefix", default = "NMask")


argv <- parse_args(p)

if(!dir.exists(argv$output)){# If the dir does not exist
  cat(paste0("Creating output directory: ", argv$output), "\n")
  dir.create(argv$output)
}


# Reads the ref
ref <- readDNAStringSet(argv$ref)
# Reads the bed 

files <- list.files(argv$bed, pattern = ".bed")

bed_list <- data.frame(names = str_remove(files, ".bed"),files)

bed_files <- list()
for(i in 1:length(files)){
  dat <- read_delim(paste0(argv$bed,"/", bed_list$files[i]), col_names = FALSE, show_col_types = FALSE)
  names(dat) <- c("chrom", "chromStart", "ChromEnd", "Name")
  dat <- mutate(dat, ref = str_remove_all(Name, "to[A-Z]"),
           alt = str_remove_all(Name, "[A-Z]to")) %>% 
    select(-Name)
  
  bed_files[[i]] <- dat
  names(bed_files)[i] <- bed_list$names[i]
}

for(i in 1:length(bed_files)){
  if(i==1){
    dat_df <- data.frame(bed_files[[1]])
    names(dat_df)[names(dat_df)=="alt"] <- bed_list$names[i]
  }else{
    dat_df <- full_join(dat_df, bed_files[[i]], by = c("chrom","chromStart","ChromEnd", "ref"))
    names(dat_df)[names(dat_df)=="alt"] <- bed_list$names[i]
  }
}

write_delim(dat_df, paste0(argv$output, "/", argv$name, "_full.bed"), delim = "\t")
dat_df2 <-dat_df %>% select(chrom, chromStart, ChromEnd)
write_delim(dat_df2, paste0(argv$output, "/", argv$name, "_compact.bed"), delim = "\t", col_names = FALSE)


system(paste0("bedtools maskfasta -fi ", argv$ref, " -bed ", paste0(argv$output, "/", argv$name, "_compact.bed"),  " -fo ", argv$output,"/",argv$name, ".fasta"))
cat(paste0("Complete \n"))

