# CURRENTLY UNDER ACTIVE DEVELOPMENT
## This README is not up to date. 

The scripts work although are unpolished, and not assembled. 
If you have an intrest in using this software, please reach out.
cxk853@student.bham.ac.uk


# NMaskGen
Generates a .fasta, were any base position with sequence variability is masked via an N

# Installation 
All R scripts needs to be made executable

chmod +x {script.R}

# Useage

Run the package installing script:

./package_install.R

To run use:

./n_mask_gen.R -s -t -n -o -m

For help use:

./n_mask_gen.R help

Example:

./n_mask_gen.R -s test_data/test_seq.fasta -t 1 -n test_output -o test_output




