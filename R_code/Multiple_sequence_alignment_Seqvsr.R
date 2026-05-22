## Author: Ana Rodríguez López
## Date: 22/05/25

# Plotting Hornet and Hormin multiple alignment with seqvisr. Use a multiple sequence alignment file in FASTA format.

# HORNET

library("seqvisr")
# Loading sample data.
testmsa <- "Unique_Hornet_genomic_seqs_Foxy3_MAFFT.fasta"

# Basic visualization with the sequence "scaffold_13_116293_119287_Hornet" set as the reference.
msavisr(mymsa = testmsa, myref = "scaffold_13_116293_119287_Hornet", basecolors = c("#000080","#D55E00","#C0C0C0"))
                                       

# HORMIN

#Loading sample data.
testmsa <- "Unique_Hormin_genomic_seqs_Foxy3_MAFFT.fasta"

#Basic visualization with the sequence "scaffold_1_5762835_5763593_Hormin" set as the reference.
msavisr(mymsa = testmsa, myref = "scaffold_1_5762835_5763593_Hormin", basecolors = c("#000080","#D55E00","#C0C0C0"))
