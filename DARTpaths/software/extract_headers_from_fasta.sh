#!/bin/bash

### this is a script to extract lines which begin with ">" character from a fasta file
### and then remove ">" from each line leaving only the entries.
### this resulting file is used as a list for all the protein IDs in that set

pathway_name=$1
protein_name=$2
### usage example below

#sed -n '/>/p' ./Microtubule_TUBA_reactome.fa > test.txt
#sed 's/[>]//g' test.txt > Microtubule_TUBA_human_protein_ID.txt


sed -n '/>/p' $HOME/DARTpaths/Proteins/${pathway_name}/${pathway_name}_${protein_name}_reactome.fa >  temp.txt
sed 's/[>]//g' temp.txt > $HOME/DARTpaths/Proteins/${pathway_name}/${pathway_name}_${protein_name}_human_protein_ID.txt

rm temp.txt
