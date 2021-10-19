#!/bin/bash

set -e

prot=$1
current_path=$(pwd)
#path to reconcile_scripts folder
path_to_folder=$2

### script to retrive ortholog from Notung output

## For Mac
#PATH=$PATH:$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/software/reconcile_scripts
### make folder for ortholog Notung_files
mkdir -p ./${prot}_orthologs
cp ./Notung_files/003_Notung_rearrange/${prot}_G_TREE_ReadyForNotung.nw.gtpruned.rooting.0.rearrange.0.homologs.txt \
./${prot}_orthologs/

echo "generating ortholog file...."
python ${path_to_folder}/005_get_orthologs.py ${prot} ${current_path}

echo "Done! see ortholog folder."
