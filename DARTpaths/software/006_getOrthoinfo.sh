#!/bin/bash

set -e

prot=$1
current_path=$(pwd)
query_protein_ID=$2

# path to reconcile_scripts folder
path_to_folder=$3

#path to folder which contains ENSEMBL databases of species
path_to_db_root=$4

### script to retrive ortholog from Notung output
### make folder for ortholog Notung_files
mkdir -p ./${prot}_orthologs/${prot}_orthoinfo

cp ./${prot}_orthologs/${prot}_human_orthologs.txt ./${prot}_orthologs/${prot}_orthoinfo/

echo "Generating orthoinfo file...."
python ${path_to_folder}/006_get_orthoinfo.py \
${prot} ${current_path} ${path_to_db_root} ${query_protein_ID}

echo "Done! check orthoinfo folder."
