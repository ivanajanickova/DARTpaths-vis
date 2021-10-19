#!/usr/local/bin/Rscript --vanilla
args <- commandArgs(trailingOnly = T)

library('ape')
library('ade4')
library('dispRity')
library('phylobase')
library('phylotools')
library('adephylo')
library('stringr')

sp_tree_001 <-  read.tree(args[1])

poly_resolved_sp_tree <- multi2di(sp_tree_001)

tree_rem_zero <- dispRity::remove.zero.brlen(poly_resolved_sp_tree, 1,verbose = FALSE)

plot(tree_rem_zero)  # this will create a pdf of poly resolved species tree in the folder
write.tree(tree_rem_zero,args[2])