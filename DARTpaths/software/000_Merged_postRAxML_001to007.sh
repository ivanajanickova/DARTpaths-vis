#!/bin/bash

set -e
current_path=$(pwd)

####################### BELOW NEEDS TO BE SPECIFIED FOR THE PROTEIN AND YOUR SYSTEM ######################################

protein_name=$1 ### e.g) AHR_ARNT  ###protein name in the format pathwayName_proteinName
## protein_ID of that protein
protein_ID=$2  #### e.g) "9606.ENSP00000351407"

###################### SPECIFY ABOVE FOR YOUR PROTEIN ###################################################################

## specify folder path where bash files are stored.
path_to_script_files=$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/software

## path to reconcile_scripts folder
path_to_reconcile_scripts=$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/software/reconcile_scripts/

## path to Notung folder.
path_to_notung=$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/software/Notung/Notung-2.9.1.3/
## if you use a different version of Notung,
## make sure to update the 004_executeNotung.sh script.
## it uses version specific notung.jar

## path to species databases folder
path_to_databases=$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/databases/

########################### THE ABOVE NEEDS TO BE SPECIFIED IN YOUR SYSTEM ################################################

## activate conda python environment which has python3.5 (named "py35" in this case)
source activate py35

echo "this is the merged script run immediately after RAxML-NG version 0.9.0 on server."
echo "running this script is the last step of the protein-phylogeny pipeline..."
echo "it will create folders for #1. orthologs #2. ortho_info files of the protein when done."
echo "this step might take a few minutes... log will be written to .log file."

# 001_postRAxML.sh is the first script to run after finishing RAxML-NG ver0.9.0 on server
bash ${path_to_script_files}/001_postRAxML.sh ${protein_name} ${path_to_reconcile_scripts} > ./${protein_name}_log_file_post-RAxML.txt
## 002_postRAxML.sh is for resolvePolytomy.R and does everything related to resolve polytomy
bash ${path_to_script_files}/002_postRAxML.sh ${protein_name} ${path_to_reconcile_scripts} >> ./${protein_name}_log_file_post-RAxML.txt
bash ${path_to_script_files}/003_postRAxML.sh ${protein_name} ${path_to_reconcile_scripts} >> ./${protein_name}_log_file_post-RAxML.txt
bash ${path_to_script_files}/004_executeNotung.sh ${protein_name} ${path_to_notung} >> ./${protein_name}_log_file_post-RAxML.txt
bash ${path_to_script_files}/005_getOrthologs.sh ${protein_name} ${path_to_reconcile_scripts} >> ./${protein_name}_log_file_post-RAxML.txt
bash ${path_to_script_files}/006_getOrthoinfo.sh ${protein_name} ${protein_ID} ${path_to_reconcile_scripts} ${path_to_databases} >> ./${protein_name}_log_file_post-RAxML.txt

bash ${path_to_script_files}/007_map_ID_for_new_orthologs.sh ${protein_name} ${protein_ID} >> ./${protein_name}_log_file_post-RAxML.txt

echo "DONE!! orthoinfo file is ready!! check folder."
echo "log written to ##protein_name##_log_file_post-RAxML.txt  exiting..."
