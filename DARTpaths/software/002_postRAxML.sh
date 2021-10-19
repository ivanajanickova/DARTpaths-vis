#!/bin/bash
set -e

prot=$1
current_path=$(pwd)
#replace_blank=$(tr -s '[:blank:]' '_' < ./${prot}_post_RAxML_processing/${prot}_polytomous_species_tree.nw)
path_to_folder=$2

source activate py35
## Assumption1: this script will be execueted
## from the /Prot/prot_family_name/protein_name/ folder
## which means, cd command should show the
## path $HOME/DARTpaths/Prot/prot_family_name/protein_name/
## this script itself must be located in
## $HOME/DARTpaths/software/ whatever system may be used.

##### script to generate species tree from list of
##### species(NCBI taxonomy ID) commandline using ete3
##### step1. this will create a polytomous species tree.
#####        this tree can be viewed with command: ete3 view -t TREE_file_NAME
##### step2. create a polytomy resolved tree (R script)

echo "creating species tree from list of NCBI taxonomy IDs..."
### point to polytomous species tree created from a list of NCBI taxID found in RAxML gene tree

### create polytomous species tree from taxonomy ID of species
cut -f1 ./${prot}_post_RAxML_processing/${prot}_species_NCBI_ID.txt | ete3 ncbiquery --tree > ./${prot}_post_RAxML_processing/${prot}_polytomous_species_tree.nw
echo "species tree ready!"
###### can be created from species name list as well
## cut -f1 ./${prot}_treespecies.txt | ete3 ncbiquery --tree > ./${prot}_post_RAxML_processing/${prot}_polytomous_species_tree.nw

#### before we give this tree file to R, we need to replace all space " " characters into something else...
### replace all the blank spaces in the file to "_" and after processing,
### that is, resolving all the polytomies, replace all the "_" to blank spaces

#### step1 ####
#### change ' ' to '_'
replace_blank=$(tr -s '[:blank:]' '_' < ./${prot}_post_RAxML_processing/${prot}_polytomous_species_tree.nw)

echo "$replace_blank" > ./${prot}_polytomous_species_tree_tmp.nw

#### step2 ####
#### create a polytomy resolved tree (R script)
chmod +rw ./${prot}_polytomous_species_tree_tmp.nw
chmod +x ${path_to_folder}/002_resolvePolytomy.R
Rscript ${path_to_folder}/002_resolvePolytomy.R ./${prot}_polytomous_species_tree_tmp.nw  ./${prot}_poly_RESOLVED_species_tree_tmp.nw

#### step3 ####
#### change '_' to space ' '
replace_underscore=$(tr -s '_' '[:blank:]' < ./${prot}_poly_RESOLVED_species_tree_tmp.nw )
echo "$replace_underscore" > ./${prot}_poly_RESOLVED_species_tree.nw

#[ADD THIS TO ANOTHER SCRIPT]##### For some species tree there are ":" and "." in the file. Remove these.

#sed 's/\./''/g' ./${prot}_poly_RESOLVED_species_tree_0.nw > ./${prot}_poly_RESOLVED_species_tree_1.nw
#sed 's/\:/''/g' ./${prot}_poly_RESOLVED_species_tree_1.nw > ./${prot}_poly_RESOLVED_species_tree.nw

echo "polytomy resolve complete!! see folder...check all species are included in the tree.. exiting..."
