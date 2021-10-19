#!/bin/bash
# -*- coding: utf-8 -*-

set -e
current_path=$(pwd)


##### USAGE #####
##### REQUIRES INTERNET CONNECTION !

###  (py35) Ds-MacBook-Air:Pathway_files d$ pwd
###  /Users/d/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/data_preparation/Pathway_files
###  bash ../00_prepare_data.sh AHR

### make sure to use the correct envirnment(python dependencies) 
##activate conda environment
#source activate py35

### specify name for pathway. use abbreviation if actual name is too long.
pathway_name=$1

## specify folder path where 00_get_uniprot_data.py file is stored. usually /DARTpaths/software
path_to_script_files=$HOME/mock_pipeline/DARTpaths/software

path_to_pathway_files=$HOME/mock_pipeline/DARTpaths

### name of file downloaded from Reactome
### keep this under /DARTpaths/Pathway_files
### can be modified to take from commandline arg
prots_file_name='AHR_pathway_proteins_R-HSA-8937144.tsv'

##### EXAMPLE OF WHAT TO SPECIFY
##  pathway_name='AHR'
##  folder_path_root='/Users/d/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/data_preparation'
##  prots_file_name='AHR_pathway_proteins_R-HSA-8937144.tsv'


python ${path_to_script_files}/00_get_uniprot_data.py \
${pathway_name} ${path_to_pathway_files} ${prots_file_name}
