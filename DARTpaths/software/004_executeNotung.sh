#!/bin/bash

set -e

prot=$1
current_path=$(pwd)
path_to_notung_sw_folder=$2

### this script is executed as:
### (py35) Ds-MacBook-Air:AHR_ARNT d$   #### <--- source activate py35  (use conda environment with python3.5)
### execute in folder of that specific protein
### /path to this script/this-script-name.sh protein_name
####(example) bash /Users/d/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/software/004_executeNotung.sh AHR_ARNT


### script to execute Notung. divided into three parts.

### step1. execute Notung's --gtpruned (prune gene tree)
### step2. execute Notung's --root (root tree)
### step3. execute Notung's --rearrange

## make directories for each of the above steps
mkdir -p ./Notung_files/001_Notung_GT_prune
mkdir -p ./Notung_files/002_Notung_root
mkdir -p ./Notung_files/003_Notung_rearrange

## copy protein tree and species tree to folder where Notung step 1 will store files
cp ./${prot}_G_TREE_ReadyForNotung.nw ./Notung_files/001_Notung_GT_prune/
cp ./${prot}_SP_TREE_ReadyForNotung.nw ./Notung_files/001_Notung_GT_prune/  ### just copying to folder...

### step 1 ###
##### --gtprune from Notung
java -jar ${path_to_notung_sw_folder}/Notung-2.9.1.3.jar \
-g ./Notung_files/001_Notung_GT_prune/${prot}_G_TREE_ReadyForNotung.nw \
-s ./${prot}_SP_TREE_ReadyForNotung.nw \
--outputdir ./Notung_files/001_Notung_GT_prune/ \
--speciestag postfix \
--edgeweights name \
--gtpruned \
--resolve \
--log

#-g ./${prot}_G_TREE_ReadyForNotung.nw \
#-s ./${prot}_SP_TREE_ReadyForNotung.nw  \

### step 2 ###
### copy pruned gene tree to folder 002_Notung_root for next step
cp ./Notung_files/001_Notung_GT_prune/${prot}_G_TREE_ReadyForNotung.nw.gtpruned ./Notung_files/002_Notung_root/
### rooting in NOTUNG ###
##### better way to specify path needed...
java -jar ${path_to_notung_sw_folder}/Notung-2.9.1.3.jar \
-g ./Notung_files/002_Notung_root/${prot}_G_TREE_ReadyForNotung.nw.gtpruned \
-s ./${prot}_SP_TREE_ReadyForNotung.nw  \
--outputdir ./Notung_files/002_Notung_root/ \
--rootscores \
--speciestag postfix \
--root \
--log

### step 3 ###
### copy rooted tree to folder 003_Notung_rearrange
cp ./Notung_files/002_Notung_root/${prot}_G_TREE_ReadyForNotung.nw.gtpruned.rooting.0 ./Notung_files/003_Notung_rearrange/

##### better way to specify path needed...
java -jar ${path_to_notung_sw_folder}/Notung-2.9.1.3.jar \
-g ./Notung_files/003_Notung_rearrange/${prot}_G_TREE_ReadyForNotung.nw.gtpruned.rooting.0 \
-s ./${prot}_SP_TREE_ReadyForNotung.nw  \
--outputdir ./Notung_files/003_Notung_rearrange/ \
--rearrange \
--treestats \
--threshold 90% \
--maxtrees 10 \
--events \
--parsable \
--savepng \
--homologtabletabs \
--log


### additional step: keep html format in designated folder ...just for convenience...
### the only difference from above rearrange is --homologtablehtml and --outputdir. same otherwise.
mkdir -p ./Notung_files/003_Notung_rearrange/html_format

java -jar ${path_to_notung_sw_folder}/Notung-2.9.1.3.jar \
-g ./Notung_files/003_Notung_rearrange/${prot}_G_TREE_ReadyForNotung.nw.gtpruned.rooting.0 \
-s ./${prot}_SP_TREE_ReadyForNotung.nw  \
--outputdir ./Notung_files/003_Notung_rearrange/html_format/ \
--rearrange \
--treestats \
--threshold 90% \
--maxtrees 10 \
--events \
--parsable \
--savepng \
--homologtablehtml \
--log
