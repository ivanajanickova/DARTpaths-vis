#!/bin/bash
set -e

#sudo apt install rename

##### following need to be specified from the command line #####
pathway_name=$1 ## e.g) AHR
protein_name=$2 ## e.g) AHR_ARNT


##### following need to be specified from the command line #####

path_to_sw_folder=$HOME/DARTpaths/software
path_to_fasta_file=$HOME/DARTpaths/Proteins/${pathway_name}/${pathway_name}_${protein_name}/${pathway_name}_${protein_name}_alignments
path_to_dhhc_list=$HOME/DARTpaths/List_of_species/Subset_species_DHHC_filtered_ID.txt
path_to_preprocess_folder=${path_to_fasta_file}/pre_RAxML_processing

## change filename sinc mv command does not accept multiple "."(extensions)
## in a filename

rename 's/\./_/g' *.*.phy
rename 's/\./_/g' *.*.log

mv -f *_DHHC_filtered_raxml_log ${path_to_fasta_file}/pre_RAxML_processing
mv -f *_DHHC_filtered_raxml_reduced_phy ${path_to_fasta_file}/pre_RAxML_processing

echo "check the log file generated in the folder..."

echo "this will tell you if the fasta file contains duplicates..."
echo "the resulting .phy file is a file with duplicates removed..."
echo "therefore it is important to check so that selected species entries are not omitted..."
echo ".reduced.phy file is the one which will be passed on to RAxML"
echo "before that, check with --parse command for compatibility.."
