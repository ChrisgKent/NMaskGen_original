#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(argparser))

# This can generate the psuedo genome from the bed file pdroduced via pseudo_genome_gen.R
## Currently being used as a sanity checker

p <- arg_parser("De novo pseudo genome generator", hide.opts = TRUE) 
# Adding flags
p <- add_argument(p, "--ref", help="Directory of the reference genome",  default = "resources/MN908947.3.fasta")
p <- add_argument(p, "--bed1", help= "Directory that contains the .bed files")
p <- add_argument(p, "--bed2", help= "(optional): Directory of a second .bed file", default = "NA")

p <- add_argument(p, "--prefix_bed1", help= "(optional): Prefix added to names from bed1", default = "NA")
p <- add_argument(p, "--prefix_bed2", help= "(optional): Prefix added to names from bed2", 
                  default = "NA")

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
files <- list.files(argv$bed1, pattern = ".bed") # Detects all files with .bed in the name
bed_list <- data.frame(names = str_remove(files, ".bed"),files) 

bed_files <- list() # Reads each file in
for(i in 1:length(files)){
  # Reads the file 
  dat <- read_delim(paste0(argv$bed1,"/", bed_list$files[i]), col_names = FALSE, show_col_types = FALSE)
  # Adds colnames
  names(dat) <- c("chrom", "chromStart", "chromEnd", "name")
  # Parses the data to generate the base change
  dat <- mutate(dat, ref = str_remove_all(name, "to[A-Z]"),
           alt = str_remove_all(name, "[A-Z]to")) %>% 
    select(-name)
  # Saves the readfile into the global Enviroment 
  bed_files[[i]] <- dat
  
  # If a prefix is given it's added to the file name 
  if(argv$prefix_bed1 == "NA"){
    names(bed_files)[i] <- bed_list$names[i]
  }else{
    names(bed_files)[i] <- paste0(argv$prefix_bed1, "_", bed_list$names[i])
    }
 
}
cat(paste0(length(files), " .bed files read in for bed1 \n"))


# If a second bed file is given
if(argv$bed2 != "NA"){
  # Finds all .bed files in the given dir
  files2 <- list.files(argv$bed2, pattern = ".bed")
  bed_list2 <- data.frame(names = str_remove(files2, ".bed"),files2)
  for(i in 1:length(files2)){
    # Reads the bed files in 
    dat <- read_delim(paste0(argv$bed2,"/", bed_list2$files[i]), col_names = FALSE, show_col_types = FALSE)
    names(dat) <- c("chrom", "chromStart", "chromEnd", "name")
    # Parses the .bed file
    dat <- mutate(dat, ref = str_remove_all(name, "to[A-Z]"),
                  alt = str_remove_all(name, "[A-Z]to")) %>% 
      select(-name)
    
    bed_files[[length(files)+i]] <- dat
    if(argv$prefix_bed2 == "NA"){
      names(bed_files)[length(files)+i] <- bed_list2$names[i]
    }else{
      names(bed_files)[length(files)+i] <- paste0(argv$prefix_bed2,"_", bed_list2$names[i])
    }
  }
  cat(paste0(length(files2), " .bed files read in for bed2 \n"))
}

# Added a warning if more than one file has the same name
if(length(unique(names(bed_files))) != length(names(bed_files))){
  stop("The same filename has been found in bed1 and bed2\nPlease rerun with something in one, or both of the --prefix_bed1 or --prefix_bed2 arguments")
}

for(i in 1:length(bed_files)){
  if(i==1){
    dat_df <- data.frame(bed_files[[1]])
    names(dat_df)[names(dat_df)=="alt"] <- names(bed_files)[i]
  }else{
    dat_df <- full_join(dat_df, bed_files[[i]], by = c("chrom","chromStart","chromEnd", "ref"))
    names(dat_df)[names(dat_df)=="alt"] <- names(bed_files)[i]
  }
}

dat_df <- arrange(dat_df, chromStart)
# Writes a file that contains the base change provided by each pseudoref
write_delim(dat_df, paste0(argv$output, "/", argv$name, "_log"), delim = "\t")

# Writes a .bed file 
dat_df2 <-dat_df %>% select(chrom, chromStart, chromEnd)
write_delim(dat_df2, paste0(argv$output, "/", argv$name, ".bed"), delim = "\t", col_names = FALSE)

# Calls bedtools maskfasta on the referance genome (-fi), provides the bed file, that contains the varied positions 
## and an output location (-fo)
system(paste0("bedtools maskfasta -fi ", argv$ref, " -bed ", paste0(argv$output, "/", argv$name, ".bed"),  " -fo ", argv$output,"/",argv$name, "_NMasked.fasta"))

# As the name of the NMasked fasta will be that of the referance genome
## The output from bedtools is read in, the name changed and resaved.
tmp <- readDNAStringSet(paste0(argv$output,"/",argv$name, "_NMasked.fasta"))
names(tmp) <- paste0(argv$name, "_NMasked")
writeXStringSet(tmp, paste0(argv$output,"/",argv$name, "_NMasked.fasta"))


cat(paste0("Compelete \n"))

