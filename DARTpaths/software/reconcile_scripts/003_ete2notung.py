#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 21:28:20 2019

@author: d
"""

#### script to prepare species tree for NOTUNG
### takes species file from resolve_polytomy.py (<--- use format 9)output
###

import sys
from Bio import Phylo
import ete3
from  ete3 import NCBITaxa, Tree, PhyloTree

protein_name=sys.argv[1]
current_path=sys.argv[2]

ncbi = NCBITaxa()

Dir_path = current_path+'/'                     # directory path root
path_to_notung_comp_species_tree_file = Dir_path+protein_name+'_SP_TREE_ReadyForNotung.nw'  # output of this script
#path_to_species_tree_file = Dir_path+protein_name+'_poly_RESOLVED_SP_tree_FRMT_9.nw'        # input to this script # 2020/03/16 commented out
path_to_species_tree_file = Dir_path+protein_name+'_poly_RESOLVED_species_tree.nw'  

poly_resolved_species_tree_file = path_to_species_tree_file
species_tree_file = Phylo.read(poly_resolved_species_tree_file,'newick')

with open(poly_resolved_species_tree_file,'r') as stree:
    text = stree.read()
#    text = text.strip()
text = text.strip()
lines = text.split(',')

rearranged_sp_list_0 = []
rearranged_sp_list = []

for n in lines:

    n_0 = ''.join(i for i in n if not i.isdigit())
    n_01 = n_0.replace('-','')


    n_1 = ''.join(word[0].upper() + word[1:] for word in n_01.split())
    rearranged_sp_list_0.append(n_1)
#print(rearranged_sp_list_0)

for sp in rearranged_sp_list_0:
    sp_1 = sp+','
    rearranged_sp_list.append(sp_1)
#print(rearranged_sp_list)

sp_str_1 = "".join(rearranged_sp_list)
#print(sp_str_1)

if sp_str_1[-2:] == ";,":
    sp_str_0 = sp_str_1.strip(',')
#print(sp_str_0)
#print(len(sp_str))


#remove last ')' and first '(' in file
# NOTUNG example species tree don't have it.

if sp_str_0[-2:] == ");":
    sp_final = sp_str_0[1:-2]+';'
#print(sp_final)
#print(sp_str_final)

with open(path_to_notung_comp_species_tree_file,"w+") as stree_ncomp:
#    stree_ncomp.write(sp_str_0)
    stree_ncomp.write(sp_final)
