# NMaskGen
Generates a .fasta, were any base position with sequence variability is masked via an N

# Installation 
Currently the dependencies need to be manually installed;

R

CRAN
-	tidyverse
-	argparser

Bioconductor
-	msa
-	biostrings

bedtools 

# Useage

You need to make the script executable. 
Move into the NMaskGen directory

chmod +x n_mask_gen.R

To run use; 
./n_mask_gen.R -s -t -n -o -m

For help use;
./n_mask_gen.R help




