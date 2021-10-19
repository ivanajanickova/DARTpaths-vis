#!/bin/bash
set -e

prot=$1
current_path=$(pwd)
path_to_folder=$2


## Assumption1: this script will be execueted
## from the /Prot/prot_family_name/protein_name/ folder
## which means, cd command should show the
## path $HOME/DARTpaths/Prot/prot_family_name/protein_name/
## this script itself must be located in
## $HOME/DARTpaths/software/ whatever system may be used.



### Step1.
## generate Notung compatible species tree from polytomy resolved species tree (tree format9 in ete3)
echo "generating Notung compatible species tree..."
python ${path_to_folder}/003_ete2notung.py ${prot} ${current_path}
echo "DONE! check for Notung compatible species tree file in folder."

###### For some species tree there are ":" and "." in the file. Remove these.

sed 's/\./''/g' ./${prot}_SP_TREE_ReadyForNotung.nw > ./${prot}_SP_TREE_ReadyForNotung_0.nw
sed 's/\:/''/g' ./${prot}_SP_TREE_ReadyForNotung_0.nw > ./${prot}_SP_TREE_ReadyForNotung.nw

### Step2.
##
echo "generating Notung compatible protein phylogeny..."
python ${path_to_folder}/003_RAxML2Notung.py ${prot} ${current_path}
echo "DONE! check for Notung compatible protein(GTREE) tree file in folder."
