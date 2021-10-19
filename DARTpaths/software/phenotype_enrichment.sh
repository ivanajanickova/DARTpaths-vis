#!/bin/bash


##### Expects to be run under pathway folder e.g.) /DARTpaths/Pathways/AHR
set -e

pathway_name=$1 ## like "Hedgehog"
Reactome_pathway_ID=$1

current_path=$(pwd)

path_to_script=$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/

#path_to_phenotype_enrichment_result=$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/Proteins/${pathway_name}/RESULT_${pathway_name}_${species}_phenotype_enrichment.txt


path_to_database_folder=$HOME/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/databases


path_to_phenotype_database_folder=${path_to_database_folder}/phenotype
path_to_ontology_database_folder=${path_to_database_folder}/ontology


echo "starting phenotype enrichment..."


python ${path_to_script}/software/phenotype_enrichment.py ${Reactome_pathway_ID} ${pathway_name} ${path_to_script}

echo "Phenotype enrichment done!! check file!"
