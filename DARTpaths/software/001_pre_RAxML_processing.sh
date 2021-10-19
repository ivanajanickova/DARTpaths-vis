#!/bin/bash
set -e

#sudo apt install rename
##### following need to be specified from the command line #####
pathway_name=$1 ## e.g) AHR
protein_name=$2 ## e.g) AHR_ARNT

path_to_sw_folder=$HOME/DARTpaths/software
path_to_fasta_file=$HOME/DARTpaths/Proteins/${pathway_name}/${pathway_name}_${protein_name}/${pathway_name}_${protein_name}_alignments
path_to_dhhc_list=$HOME/DARTpaths/List_of_species/Subset_species_DHHC_filtered_ID.txt
path_to_preprocess_folder=${path_to_fasta_file}/pre_RAxML_processing

echo "starting step1"
echo "first it will filter the FASTA file(the one which ends with .aln-gb in ##protein-name##_alignments folder ) "
echo "so that it only includes data for the species included in DHHC list of species."

#### STEP1 ####
mkdir -p ${path_to_fasta_file}/pre_RAxML_processing

## create a blank file that will store matching headers for this protein.
touch ${path_to_fasta_file}/pre_RAxML_processing/${pathway_name}_${protein_name}_matching_headers_tmp.txt
touch ${path_to_fasta_file}/pre_RAxML_processing/${pathway_name}_${protein_name}_matching_headers.txt
touch ${path_to_fasta_file}/pre_RAxML_processing/${pathway_name}_${protein_name}fam_mafft_hmm_DHHC_filtered.aln-gb

## extract items from fasta file(which ends with .aln-gb This.aln-gb file is an output from reconcile1.sh)
## species which are in DHHC species list.
## with this, the fasta file will contain data only for species in DHHC list  (done in next step)

fgrep -f \
${path_to_dhhc_list} \
${path_to_fasta_file}/${pathway_name}_${protein_name}fam_mafft_hmm.aln-gb > ${path_to_preprocess_folder}/${pathway_name}_${protein_name}_matching_headers_tmp.txt

#### remove ">" character from the begining

sed 's/[>]//g' ${path_to_preprocess_folder}/${pathway_name}_${protein_name}_matching_headers_tmp.txt > ${path_to_preprocess_folder}/${pathway_name}_${protein_name}_matching_headers.txt

${path_to_sw_folder}/faSomeRecords ${path_to_fasta_file}/${pathway_name}_${protein_name}fam_mafft_hmm.aln-gb ${path_to_preprocess_folder}/${pathway_name}_${protein_name}_matching_headers.txt \
 ${path_to_preprocess_folder}/${pathway_name}_${protein_name}fam_mafft_hmm_DHHC_filtered.aln-gb

##### STEP2 #####
echo "done!! see folder for ##protein-name##fam_mafft_hmm_DHHC_filtered.aln-gb file "
echo "this will be passed on to process for input to RAxML-NG (version 0.9.0)"
echo "starting step2. preprocessing for RAxML-NG(version 0.9.0)"

#cd ${path_to_sw_folder}/RAxML-NG_0-9-0
echo "first it will create a temp file  "
${path_to_sw_folder}/RAxML-NG_0-9-0/raxml-ng --check --msa ${path_to_preprocess_folder}/${pathway_name}_${protein_name}fam_mafft_hmm_DHHC_filtered.aln-gb --model LG+FC+G8m --prefix ${protein_name}_DHHC_filtered

## RAxML-NG creates two files #1) .phy file and #2) log file in the current folder.
## move these to our working folder. done in next script 002_pre_RAxML_processing.sh
