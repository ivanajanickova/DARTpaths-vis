#!/bin/bash
set -e

prot=$1
## path to reconcile_scripts folder. this script is called from 000_Merged_001to006.sh
path_to_folder=$2
current_path=$(pwd)

source activate py35

## Assumption1: this script will be execueted
## from the /Prot/prot_family_name/protein_name/ folder
## which means, cd command should show the
## path $HOME/DARTpaths/Prot/prot_family_name/protein_name/
## this script itself must be located in
## $HOME/DARTpaths/software/ whatever system may be used.

## For Mac
#PATH=$PATH:$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/software/reconcile_scripts

echo "This is the first script to be run after RAxML-NG version 0.9.0. It will create two files:"
echo "#1.a list of species in the protein tree, #2.polytomous species tree file."

mkdir -p ${prot}_post_RAxML_processing

## copy the protein phylogeny file which ends with
## .support (tree file with bootstrap support values) to
##  the above folder. ${prot}_post_RAxML_processing
cp ./${prot}_trees/*.support ./${prot}_post_RAxML_processing/

### to do: instead of specifying the entire path to the python script file like below, retrieve from PATH above.

python ${path_to_folder}/001_for_raxml_ver9_specnr2specnames.py ${prot} ${current_path}
## python 001_for_raxml_ver9_specnr2specnames.py ${prot} ${current_path}

echo "DONE! check folder for file of species list that are in the protein phylogeny..."

taxonkit name2taxid ./${prot}_post_RAxML_processing/${prot}_treespecies.txt > ./${prot}_post_RAxML_processing/${prot}_treespecies_name_and_NCBI_ID.txt

##[NOT IN USE: gives error]grep -Eo '[0-9\.]+' ./${prot}_post_RAxML_processing/${prot}_treespecies_name_and_NCBI_ID.txt > ./${prot}_post_RAxML_processing/${prot}_species_NCBI_ID.txt

### From ${prot}_treespecies_name_and_NCBI_ID.txt file, extract the last column and write to file ${prot}_post_RAxML_processing/${prot}_species_NCBI_ID.txt
### this assumes that the last column contains (which is usually the case) "species_NCBI_ID"

awk '{print $NF}' ./${prot}_post_RAxML_processing/${prot}_treespecies_name_and_NCBI_ID.txt > ./${prot}_post_RAxML_processing/${prot}_species_NCBI_ID_tmp.txt

### some species have non-numerical values which look like (e.g.) "CD1", "DSM44101". our list of eukaryote species doesn't contain these.
### Therefore, we will exclude any line which has non-integers. Which leaves only the relevant ones. (avoid error)

grep -E '^[[:digit:]]+$' ./${prot}_post_RAxML_processing/${prot}_species_NCBI_ID_tmp.txt > ./${prot}_post_RAxML_processing/${prot}_species_NCBI_ID.txt
