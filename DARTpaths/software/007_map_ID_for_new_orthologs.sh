#!/bin/bash

set -e

prot=$1   ##e.g)  "Hedgehog_UB"
#current_path=$(pwd)
query_protein_ID=$2  ### "e.g.) "9606.ENSP00000441543"

current_path=$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths

# specify where this 007_map_ID_for_new_orthologs.py script is.
software_path=${current_path}/software/reconcile_scripts

# path to pathways folder
path_to_pathway_folder=${current_path}/Proteins

## path to databases folder
path_to_database_folder=${current_path}/databases/for_id_mapping


python ${software_path}/007_map_ID_for_new_orthologs.py \
${prot} ${query_protein_ID} ${path_to_pathway_folder} ${path_to_database_folder}

echo "DONE!! check files. NOTE:look into log file. Make sure there are no errors! "
