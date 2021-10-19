#!/usr/bin/env python3

import sys

protein_name=sys.argv[1]
current_path=sys.argv[2]

# Script to get the orthologs for all human genes from the homolog table from
# Notung (Notung plays it safe and calls all them homologs)

# Input: -current dir, -name of protein, -stree number, -btree number,
#  -"all" (search all species) or leave empty to search for
# species specified in specieslist

## WARNING: input 'all' will take a considarable longer time,
#  it is almost always best to specifiy which species to look for

specieslist = ['MusMusculus','RattusNorvegicus','OryctolagusCuniculus','DanioRerio','CaenorhabditisElegans','DictoysteliumDiscoideum']

#path = sys.argv[1]
#prot = sys.argv[2]
#stree = sys.argv[3]
#btree = sys.argv[4]
#if 	len(sys.argv) > 5:
#	doit = sys.argv[5]
#else:
#doit = 'not_all'
doit = 'not_all'

## specify path to file that was output by rearrange suffix is, .homologs.txt
#path_to_notung_homologfile = "/home/d/DARTpaths/Proteins/AHR/default_Gblock_overlap_parameters/AHR2_AHRR/notung_complete_files/0930_GENE_TREE_AHR2_AHRR_for_Strip_Species.nw.gtpruned.rearrange.0.homologs.txt"
path_to_notung_homologfile = current_path+'/'+protein_name+'_orthologs'+'/'+protein_name+"_G_TREE_ReadyForNotung.nw.gtpruned.rooting.0.rearrange.0.homologs.txt"
path_to_output_file = current_path+'/'+protein_name+'_orthologs'+'/'+protein_name+"_human_orthologs.txt"


#get table, strip it from the complemantary index at top
with open(path_to_notung_homologfile,'r') as table:
	ttext = table.read()
ttext = ttext.strip()
pieces = ttext.split('\n\n')
#print(len(pieces))
#print(pieces[0])      # the first row is description
ttext = pieces[1]
#print(len(pieces[1]))

# get all genes from the first line and get a list of the genes and their P/O's
lines = ttext.split('\n')
#print(lines[1])
#print(len(lines))
#print(lines[0])

gene_index = lines[0].split('\t')
#print(len(gene_index))
#print(gene_index)

genes = lines[1:]
#print(lines[1])

human_orthos = {}

#loop over genes, if human gene: get orthologs using the gene-index,
# put everything in new list (yay more lists), combine list into string
# and add string to main-list
for gene in genes:
	tabs = gene.split('\t')
	gname = tabs[0]
#print(gname)
	if "Homo" not in gname:	# No human gene
		continue
	orthos = []
	for index, tab in enumerate(tabs):
		if tab == 'O': 	#found ortholog, so look it up and add it to list
			ortho = gene_index[index]
			pieces = ortho.rsplit('_',1)
			species = pieces[-1]
#            print(species)
			if (species in specieslist) or (doit == 'all'):
				orthos.append(ortho)
	#Now we have the orthologs, add it to the dictionary with human gene as key
	new_gene = {gname:sorted(orthos)}
	human_orthos.update(new_gene)

print('Found {} human proteins with orthologs.'.format(len(human_orthos)))
#print(human_orthos)

#save dictionary as text
with open(path_to_output_file,"w+") as save_file:
	for gene, orthos in sorted(human_orthos.items()):
		save_file.write('{}\t{}\n\n'.format(gene,', '.join(orthos)))
