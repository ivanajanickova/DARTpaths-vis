#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 23:43:03 2019

@author: d
"""

import re
import sys
import ete3
from Bio import Phylo
from  ete3 import NCBITaxa, Tree, PhyloTree

protein_name=sys.argv[1]
current_path=sys.argv[2]

#######  script to convert output from RAxML-NG version 0.9.0 's .support file for input to NOTUNG

### (***** NOT APPLICABLE *****)STEP1. convert .support to NHX format in the below link

######  http://phylogeny.lirmm.fr/phylo_cgi/get_result.cgi

### STEP2. species name from .support file is in format like  __Species_name_:

### convert this to __SpeciesName:

### STEP3. check if it is acceptable by NOTUNG.

Dir_path = current_path+'/'
gene_tree_name = Dir_path+protein_name+'_post_RAxML_processing/'+protein_name+'fam_tree.raxml.support'
gene_tree_ready_for_NOTUNG = Dir_path+protein_name+'_G_TREE_ReadyForNotung.nw'

with open(gene_tree_name,'r') as gtree:
	gtext = gtree.read()
gtext = gtext.strip()
lines = gtext.split(',')

gtree_in = Phylo.read(gene_tree_name,'newick')
#print(gtree_in)


no_of_term = gtree_in.count_terminals()  # number of terminals in gene tree from RAxML-NG
#print(no_of_term)

punc_filter = re.compile('([_]\s*)')
item_list_0=[]

####  USED FOR IDENTIFYING EDGE LENGTHS  #####

numeric_const_pattern = '(\: \d* \. \d+)'  # find all edge distances , floating point number here
rx = re.compile(numeric_const_pattern, re.VERBOSE)

for n in lines:
    n_1 = ''.join(word[0].upper() + word[1:] for word in n.split())

    first_half,last_half = n_1.split("_:")
    sp_id,sp_name_pre = first_half.split("__",1)

#    punc_filter = re.compile('([_]\s*)')
    split_with_punctuation = punc_filter.split(sp_name_pre)
    sp_name_f = ''.join([i.capitalize() for i in split_with_punctuation])
    sp_name_final = sp_name_f.replace("_","")

    edge_lengths_list = re.findall("\d+\.\d+", last_half)


    n_4 = re.sub('|'.join(edge_lengths_list), ':', last_half)
#    print(n_4)
    n_5 = re.sub(':', '',n_4)
#    print(n_5)

#
#    for i in edge_lengths_list:
#        re.sub(i,"")

#    processed_last_half = re.sub("",re.findall("\d+\.\d+", last_half)[0])

##### this n_2 only necessary if species_id needs to be re-formatted. ignore otherwise

    n_2 = sp_id.replace(".","_")

    # replace sp_id with n_2 if necessary for next steps in processing

    final_item_0 = sp_id + "__" + sp_name_final + ":" + n_5 + ","
    final_item_1 = re.sub(":\)",")",final_item_0)
    final_item = re.sub(":,",",",final_item_1)

    item_list_0.append(final_item)

#print(n_1)
#print(first_half)
#print(last_half)
#print(sp_id)
#print(sp_name_pre)
#print(split_with_punctuation)
#print(sp_name_f)
#print(sp_name_final)

#print(edge_lengths_list)

#print(type(edge_lengths_list))

#print(n_2)
#print(final_item_0)
#print(final_item)

sorted_string_0 = "".join(item_list_0)
sorted_string_1 = sorted_string_0.rstrip(',')
sorted_string_2 = sorted_string_1[1:]
#print(sorted_string)



### if error occurs in opening this file with NOTUNG, check if the parentheses are balanced.
### if Not, balance and see if it fixes.

##### check parentheses count and balance if not balanced.
### from : https://codereview.stackexchange.com/questions/153078/balanced-parentheses-checker-in-python

def add_vectors(a, b):
    for i, x in enumerate(b):
        a[i] += x  # a is a list, so it supports item assignment

    return a

def is_balanced(string):
    #         (  [  {
    counts = [0, 0, 0]

    # dict of tuples to add based on char (like a switch statement)
    d = {'(': ( 1, 0, 0), '[': ( 0, 1, 0), '{': ( 0, 0, 1),
         ')': (-1, 0, 0), ']': ( 0,-1, 0), '}': ( 0, 0,-1)}

    for char in string:
        try:
            counts = add_vectors(counts, d[char])
        except KeyError:  # char is not in dict, so it isn't '(', '{', etc.
            continue

        if -1 in counts:  # close without open
            return False

    return counts == [0, 0, 0]  # will resolve to initial state if correct

balance_check = is_balanced(sorted_string_2)
#print(balance_check)

if balance_check == False:
	sorted_string = (sorted_string_2[:-2]+';')
	### assuming that the last ')' in string is the one extra parenthesis, remove this from the right tip
	### script can be improved to be much robust....
#print(sorted_string)
##### specify  " --edgeweights name "  in first step of NOTUNG command line when using
#####   --gtpruned OR ( -resolve to convert non -binary gene tree to binary )

with open(gene_tree_ready_for_NOTUNG,"w+") as for_notung:
	for_notung.write(sorted_string)
