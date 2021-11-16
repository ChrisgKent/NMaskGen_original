#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(argparser))

# This can generate the psuedo genome from the bed file produced via pseudo_genome_gen.R
## Currently being used as a sanity checker

p <- arg_parser("De novo pseudo genome generator") 
# Adding flags
p <- add_argument(p, "--ref", help="Directory of the reference genome",  default = "resources/MN908947.3.fasta")
p <- add_argument(p, "--bed", help= "Directory of the .bed file")
p <- add_argument(p, "--output", help="Directory the de novo pseudo_genome is written", default = "FALSE")
p <- add_argument(p, "--check", help="Directory the original pseudo_genome to check against", default = "FALSE")

argv <- parse_args(p)



# Reads the ref
ref <- readDNAStringSet(argv$ref)
# Reads the bed 
bed <- read_delim(argv$bed, col_names = FALSE, show_col_types = FALSE) 
names(bed) <- c("chrom", "chromStart", "ChromEnd", "Name")

# Parses the bed file 
par_bed <- bed %>% mutate(ref = str_remove_all(Name, "to[A-Z]"),
                          alt = str_remove_all(Name, "[A-Z]to"))

# Alters the ref genome
de_novo_pseudo <- ref
names(de_novo_pseudo) <- "de_novo_pseudo"
for(i in 1:nrow(par_bed)){
  de_novo_pseudo[[1]][par_bed$ChromEnd[i]] <- par_bed$alt[i]
}

# Checks if the de novo genome is the same as the pipeline one
if(argv$check != "FALSE"){
  pipline <- readDNAStringSet(argv$check)
  if(pipline == de_novo_pseudo){cat("\nDe novo genome is identical\n")
  }else{cat("ERROR: non-identical")}
}

# Writes the de novo genome
if(argv$output != "FALSE"){
  writeXStringSet(de_novo_pseudo, argv$output)
}



