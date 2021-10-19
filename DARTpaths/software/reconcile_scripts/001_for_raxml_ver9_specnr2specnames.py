#!/usr/bin/env python3
### executed from 001_postRAxML.sh

import sys

protein_name=sys.argv[1]
current_path=sys.argv[2]

### the current_path must be set to the current directory.
### in this case, ${Pathway}/${Pathway}_${proteinname}
### e.g) ${HOME}/DARTpaths/Proteins/Angiogenesis/Angiogenesis_VEGFA

# Script takes RAxML-tree and prepare it for Notung.
# Places speciesnames in front of protein and removes species-number.
# Removes tree distances and puts bootstrap value in its place, adjusts format.

### RAxML-NG ver0.9.0 's output file (which ends with .support)

#RAxML_NG_outputfile
RAxML_NG_outputfile = current_path+'/'+protein_name+'_post_RAxML_processing'+'/'+protein_name+'fam_tree.raxml.support'

#polytomous_species_tree_file = dir_path+prot_name+'/'+sp_tree_file
polytomous_species_tree_file = current_path+'/'+protein_name+'_post_RAxML_processing'+'/'+protein_name+'_polytomous_species_tree.nw'

#e.g.) list_of_treespecies_file = '/home/d/DARTpaths/Proteins/AHR/test03_AHR5/test03_AHR5_trees/test03_AHR5fam_treespecies'
#list_of_treespecies_file = dir_path+prot_name+'/'+sp_list_file
list_of_treespecies_file = current_path+'/'+protein_name+'_post_RAxML_processing'+'/'+protein_name+'_treespecies.txt'

with open(RAxML_NG_outputfile,'r') as treefile:

    tree = treefile.read()
    tree = tree.replace('lcl|','')
    tree = tree.replace('(','(\n')
    tree = tree.replace(' (\n',' (')
    tree = tree.replace(')',')\n')
    tree = tree.replace(',',',\n')
    tree = tree.strip()

## print(tree)

# Go over lines, depending on line 'type', reformat and add line to new list, if species found add to species list
spec_list = []
newlines = []
for line in tree.split('\n'):
    if line[0] == '(':
        newlines.append(line)


    elif (line[0].isdigit()) and ('__' in line):

#        [code,name] = line.split('(')
        temp_0 = line.split('__')
        code = temp_0[0]
        temp_2 = temp_0[1]
        temp_3 = temp_2.split('_:')
        name_0 = temp_3[0]
        name = name_0.replace('_', ' ')


        if name not in spec_list:
            spec_list.append(name)
            fpoint = code.find('.')   #first occurence
            protein = code[fpoint+1:]
            newline = '{}.{}'.format(name,protein)
            newlines.append(newline)

    elif (line[0].isdigit()) and (':' in line):
##        [bs,dis] = line.split(':')

         temp_1=line.split(':')
         bs=temp_1[0]
         dis=temp_1[1]

         newlines.append(':{}{}'.format(bs,line[-1]))
    else:
        newlines.append(line[-1])

# Join new list, remove some lines, remove spaces
#newtree = '\n'.join(newlines)
#newtree = newtree.replace('\n,',',')
#newtree = newtree.replace('\n)',')')
#newtree = newtree.replace('\n;',';')
#newtree = newtree.replace(' ','_')

# Save new tree:
#with open(polytomous_species_tree_file,'w') as treefile:
#    treefile.write(newtree)

# Save species:
with open(list_of_treespecies_file,'w') as specs:
    specs.write('\n'.join(spec_list))
