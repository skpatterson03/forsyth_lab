# This file is a template for how to create motif stack plots.
# The data I am using is a csv file that only contains one type of 
# motif.  In this case I am looking at RCGDAD.  The package is ggseqlogo.  

# For more info about using this package: https://omarwagih.github.io/ggseqlogo/

# You may have to install the package on your computer

#install.packages("ggseqlogo")
library(ggseqlogo)
library(ggplot2)

# I created this file using the 'context' section in the gff
filepath = '../data/rcgdad_variants.csv'

# read in the list of RCGDAD motifs
motifs = read.csv(filepath)

# this has nice colors
ggseqlogo(motifs, method="prob", seq_type='dna')

# this is another way to produce the plot
## note: if you omit method="prob", the y-axis is in "bits"
ggplot() + geom_logo( motifs ) + theme_logo()

# this shows strong vs weak bonds
ggseqlogo(motifs, method="prob", col_scheme='base_pairing')
