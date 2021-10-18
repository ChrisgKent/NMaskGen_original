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


./n_mask_gen.R -s test_data/test_seq.fasta -t 1 -n test_output -o test_output




