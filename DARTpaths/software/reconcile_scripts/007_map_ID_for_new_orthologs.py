#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 12:46:02 2020

@author: d
"""

## TASK ### add orthologs identified by orthology mapping pipeline to phenotype enrichment set.
##### since ortholog databases only include genes/proteins that are part of ENSEMBL compara orthologs,
##### orthologs newly identified by our pipeline cannot retrieve gene_ID from this.
#####
##### Therefore, Protein_ID ---> Gene_ID mapping uses data from BioMart.

### File example:
### ENSEMBL_Mouse_Gene_Prot_ID_mapping.txt

##### First step is to list all orthologs that were identified.
##### list Human STRING_ID, Human_GENE_ID, Orhtolog_protein_ID, Ortholog_Gene_ID [per species]


### From _searches_literature folder, get STRING_ID and ENSEMBL_GENE_ID
### From ortholog result file, extract result for this STRING ID.
### divide this into per-species and from the mapping_database,
### retrieve ENSEMBL_Gene_ID for these species' Protein_ID.

### Add this into a data frame (pandas)

### [IN A SEPARATE SCRIPT] merge this with ENSEMBL_COMPARA retrieved dataframe.

### resolve duplicates and only keep unique entries if necessary.
### on this merged set, run phenotype enrichment.

import re
import os
import numpy as np
import csv
import sys
import pandas as pd
import math
import glob

###########################################################################

protein_name = sys.argv[1] ## e.g.) 'Hedgehog_UB'
query_protein_id = sys.argv[2]  ## e.g.) "9606.ENSP00000441543"
path_to_pathways_folder = sys.argv[3]


pathway_name , _ = protein_name.split('_',1) ## e.g) 'Hedgehog'

path_to_db_root = sys.argv[4]

###########################################################################

path_to_prot_folder= path_to_pathways_folder + '/'+ pathway_name+'/'+ protein_name

_ , query_protein_id_code = query_protein_id.split('.',1)
### e.g) query_protein_id_code = "ENSP00000441543"

#####
path_to_orthofiles = path_to_pathways_folder+"/"+pathway_name+"/"

#print(path_to_orthofiles )
####
#gene_ID_of_this_human_protein = "ENSG00000150991"

#human_IDs = pd.DataFrame({'Human_protein_ID':[query_protein_id_code] , 'Human_Gene_ID':[gene_ID_of_this_human_protein] })

#path_to_orthologfile = "/Users/d/Documents/PythonScripts/20191208/AHR_pathway_orthologs/"
full_path_to_ortholog_file = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+protein_name+"_orthoinfo"+'/'+protein_name+"_human_orthologs.txt"
path_to_ortho_info_file = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+protein_name+"_orthoinfo"+'/'+protein_name+"_QueryProt_orthoinfo.txt"


##### Files which contain detailed information for orthologs to query prot(matched to STRING ID 9606.ENSP ---)
##### This happens when ENSEMBL protein ID are "Retired"
#### example: http://www.ensembl.org/Danio_rerio/Transcript/Idhistory/Protein?t=ENSDARP00000089879

path_to_Mouse_info_file = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Mouse_ortho_info_"+protein_name+".csv"
path_to_Zebraf_info_file = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Zebrafish_ortho_info_"+protein_name+".csv"
path_to_Cele_info_file = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Cele_ortho_info_"+protein_name+".csv"
path_to_Rat_info_file = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Rat_ortho_info_"+protein_name+".csv"
path_to_Dicty_info_file = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Dicty_ortho_info_"+protein_name+".csv"
path_to_Rabbit_info_file = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Rabbit_ortho_info_"+protein_name+".csv"

##### Files which contain ortholog Prot ID for which Gene ID could not be found. Search on the web.
Zebraf_path_to_not_found_prot_ID = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Zebrafish_ID_NOT_FOUND_"+protein_name+".csv"
Mouse_path_to_not_found_prot_ID = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Mouse_ID_NOT_FOUND_"+protein_name+".csv"
Cele_path_to_not_found_prot_ID = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Cele_ID_NOT_FOUND_"+protein_name+".csv"
Rat_path_to_not_found_prot_ID = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Rat_ID_NOT_FOUND_"+protein_name+".csv"
Dicty_path_to_not_found_prot_ID = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Dicty_ID_NOT_FOUND_"+protein_name+".csv"
Rabbit_path_to_not_found_prot_ID = path_to_prot_folder+'/'+protein_name+"_orthologs"+'/'+"Rabbit_ID_NOT_FOUND_"+protein_name+".csv"


##### Files with ID mapping of orthologs found for this set.

Zebraf_file_with_ortho_ID = path_to_orthofiles +"Orthologs_"+pathway_name+"_pathway/"+"Zebrafish_orthologs/"+protein_name+"_orthologs.csv"
Mouse_file_with_ortho_ID = path_to_orthofiles +"Orthologs_"+pathway_name+"_pathway/"+"Mouse_orthologs/"+protein_name+"_orthologs.csv"
Cele_file_with_ortho_ID = path_to_orthofiles+"Orthologs_"+pathway_name+"_pathway/"+"Cele_orthologs/"+protein_name+"_orthologs.csv"
Rat_file_with_ortho_ID = path_to_orthofiles +"Orthologs_"+pathway_name+"_pathway/"+"Rat_orthologs/"+protein_name+"_orthologs.csv"
Dicty_file_with_ortho_ID = path_to_orthofiles +"Orthologs_"+pathway_name+"_pathway/"+"Dicty_orthologs/"+protein_name+"_orthologs.csv"
Rabbit_file_with_ortho_ID =  path_to_orthofiles +"Orthologs_"+pathway_name+"_pathway/"+"Rabbit_orthologs/"+protein_name+"_orthologs.csv"


##########

cele_path_to_db = path_to_db_root + "/ENSEMBL_Celegans_Gene_Prot_ID_mapping.txt"
rabbit_path_to_db = path_to_db_root + "/ENSEMBL_Rabbit_Gene_Prot_ID_mapping.txt"
danio_path_to_db = path_to_db_root + "/ENSEMBL_Zebrafish_Gene_Prot_ID_mapping.txt"
rat_path_to_db = path_to_db_root + "/ENSEMBL_Rat_Gene_Prot_ID_mapping.txt"
mouse_path_to_db = path_to_db_root + "/ENSEMBL_Mouse_Gene_Prot_ID_mapping.txt"
dicty_path_to_db = path_to_db_root + "/ENSEMBL_SlimeMould_Gene_Prot_ID_mapping.txt"

df_Cele = pd.read_csv(cele_path_to_db,delimiter='\t', header = 'infer')
df_Rabbit = pd.read_csv(rabbit_path_to_db,delimiter='\t', header = 'infer')
df_Danio = pd.read_csv(danio_path_to_db,delimiter='\t', header = 'infer')
df_Rat = pd.read_csv(rat_path_to_db,delimiter='\t', header = 'infer')
df_Mouse = pd.read_csv(mouse_path_to_db,delimiter='\t', header = 'infer')
df_Dicty = pd.read_csv(dicty_path_to_db,delimiter='\t', header = 'infer')

#print(df_Mouse)



with open(full_path_to_ortholog_file,'r') as ortho_tab:
	ortho_text = ortho_tab.read()

ortho_text.strip()
human_prots = ortho_text.split("\n\n")

list_of_sp_orthologs = []
list_of_NON_query_protein = []

sp_ortholog_list = []

for i in human_prots:
    if i == "":
        continue
    else:
        human_prot_full, right_side = i.split('\t')
        human_prot_id, h_sp_name = human_prot_full.split('__')

        if human_prot_id == query_protein_id:
            list_of_sp_orthologs.append(right_side)
        else:
            list_of_NON_query_protein.append(i)

#print(list_of_sp_orthologs)
#print(len(list_of_sp_orthologs))
#print(type(list_of_sp_orthologs))

sp_ortholog_list_0 = []
sp_ortholog_list = []

for i in list_of_sp_orthologs:
    i.strip()
    sp_ortholog_list_0 = i.split(',')

#print(sp_ortholog_list_0)
#print(len(sp_ortholog_list_0))
#
for i in sp_ortholog_list_0:
    i = i.strip()
    sp_ortholog_list.append(i)
#print(sp_ortholog_list)


Mouse_ortholog_prot = []
Zebrafish_ortholog_prot = []
Rat_ortholog_prot = []
Cele_ortholog_prot = []
Rabbit_ortholog_prot = []
Dicty_ortholog_prot = []

other_prot =[]

for i in sp_ortholog_list:

    if i == "":
        continue
    else:
        ncbi_id , ensmbl_prot_id_sp = i.split('.',1) #safe for cele. split at first(left-most '.')
        ensmbl_prot_id , sp_name = ensmbl_prot_id_sp.split('__',1)


        if ncbi_id == '10090':
            Mouse_ortholog_prot.append(ensmbl_prot_id)

        if ncbi_id == '6239':
            Cele_ortholog_prot.append(ensmbl_prot_id)

        if ncbi_id == '10116':
            Rat_ortholog_prot.append(ensmbl_prot_id)

        if ncbi_id == '9986':
            Rabbit_ortholog_prot.append(ensmbl_prot_id)

        if ncbi_id == '44689':
            Dicty_ortholog_prot.append(ensmbl_prot_id)

        if ncbi_id == '7955':
            Zebrafish_ortholog_prot.append(ensmbl_prot_id)

        else:
            other_prot.append(ensmbl_prot_id)

#print(Cele_ortholog_prot)
#print(Mouse_ortholog_prot)
#print(Zebrafish_ortholog_prot)
#print(Rat_ortholog_prot)

################

##### Zebrafish ######
z_not_found =[]
z_found_ortho_id = pd.DataFrame()
z_ortho_id_0 = pd.DataFrame()  ## For pheno_mapping

for i in Zebrafish_ortholog_prot:

    search_Zebraf_db = df_Danio.loc[(df_Danio['Protein stable ID'] == i) ]
#    id_list_zebrafish = search_Zebraf_db['Gene stable ID','Protein stable ID']

    if search_Zebraf_db.empty:
        z_not_found.append(i)
#        print(" The following protein cannot be found in Database. Might be an old ID. sometimes there are -Retired ENSEMBL IDS-. search manually.Record entered in NOT FOUND file. \n\n"+i)
    else:
        z_found_ortho_id = z_found_ortho_id.append(search_Zebraf_db)

z_found_ortho_id.to_csv(path_to_Zebraf_info_file)
z_ortho_id_0 = z_found_ortho_id.reindex(columns = ['Gene stable ID','Protein stable ID','Gene name'])
z_ortho_id_1 = z_ortho_id_0.rename(columns = {'Gene stable ID': 'Zebrafish_Gene_stable_ID', 'Protein stable ID':'Zebrafish_Protein_stable_ID', 'Gene name':'Zebrafish_Gene_name'})
z_ortho_ids = z_ortho_id_1.drop_duplicates() ## unique IDs of orthologous proteins

z_ortho_ids.to_csv(Zebraf_file_with_ortho_ID, header=None)

#print(z_ortho_ids)
#print(z_not_found)
#print(z_found_ortho_id)
with open(Zebraf_path_to_not_found_prot_ID, "w+") as not_found:
    wr = csv.writer(not_found)
    wr.writerow(z_not_found)

############ Zebrafish END ############


##### Mouse ######

m_not_found =[]
m_found_ortho_id = pd.DataFrame()
m_ortho_id_0 = pd.DataFrame()  ## For pheno_mapping

for i in Mouse_ortholog_prot:

    search_Mouse_db = df_Mouse.loc[(df_Mouse['Protein stable ID'] == i) ]

    if search_Mouse_db.empty:
        m_not_found.append(i)
#        print(" The following protein cannot be found in Database. Might be an old ID. sometimes there are -Retired ENSEMBL IDS-. search manually.Record entered in NOT FOUND file. \n\n"+i)
    else:
        m_found_ortho_id = m_found_ortho_id.append(search_Mouse_db)

m_found_ortho_id.to_csv(path_to_Mouse_info_file)
m_ortho_id_0 = m_found_ortho_id.reindex(columns = ['Gene stable ID','Protein stable ID','Gene name'])
m_ortho_id_1 = m_ortho_id_0.rename(columns = {'Gene stable ID': 'Mouse_Gene_stable_ID', 'Protein stable ID':'Mouse_Protein_stable_ID', 'Gene name':'Mouse_Gene_name'})
m_ortho_ids = m_ortho_id_1.drop_duplicates() ## unique IDs of orthologous proteins

m_ortho_ids.to_csv(Mouse_file_with_ortho_ID, header=None)

#print(m_ortho_ids)


with open(Mouse_path_to_not_found_prot_ID, "w+") as not_found:
    wr = csv.writer(not_found)
    wr.writerow(m_not_found)

############ Mouse END ############


###### Cele ######
#
c_not_found =[]
c_found_ortho_id = pd.DataFrame()
c_ortho_id_0 = pd.DataFrame()  ## For pheno_mapping

for i in Cele_ortholog_prot:

    search_Cele_db = df_Cele.loc[(df_Cele['Protein stable ID'] == i) ]

    if search_Cele_db.empty:
        c_not_found.append(i)
#        print(" The following protein cannot be found in Database. Might be an old ID. sometimes there are -Retired ENSEMBL IDS-. search manually.Record entered in NOT FOUND file. \n\n"+i)
    else:
        c_found_ortho_id = c_found_ortho_id.append(search_Cele_db)

c_found_ortho_id.to_csv(path_to_Cele_info_file)
c_ortho_id_0 = c_found_ortho_id.reindex(columns = ['Gene stable ID','Protein stable ID','Gene name'])
c_ortho_id_1 = c_ortho_id_0.rename(columns = {'Gene stable ID': 'Celegans_Gene_stable_ID', 'Protein stable ID':'Celegans_Protein_stable_ID', 'Gene name':'Celegans_Gene_name'})
c_ortho_ids = c_ortho_id_1.drop_duplicates() ## unique IDs of orthologous proteins

c_ortho_ids.to_csv(Cele_file_with_ortho_ID, header=None)

#print(c_ortho_ids)

with open(Cele_path_to_not_found_prot_ID, "w+") as not_found:
    wr = csv.writer(not_found)
    wr.writerow(c_not_found)
#
############# Cele END ############
#
###### Rabbit ######
#
rabbit_not_found =[]
rabbit_found_ortho_id = pd.DataFrame()
rabbit_ortho_id_0 = pd.DataFrame()  ## For pheno_mapping

for i in Rabbit_ortholog_prot:

    search_Rabbit_db = df_Rabbit.loc[(df_Rabbit['Protein stable ID'] == i) ]

    if search_Rabbit_db.empty:
        rabbit_not_found.append(i)
#        print(" The following protein cannot be found in Database. Might be an old ID. sometimes there are -Retired ENSEMBL IDS-. search manually.Record entered in NOT FOUND file. \n\n"+i)
    else:
        rabbit_found_ortho_id = rabbit_found_ortho_id.append(search_Rabbit_db)

rabbit_found_ortho_id.to_csv(path_to_Rabbit_info_file)
rabbit_ortho_id_0 = rabbit_found_ortho_id.reindex(columns = ['Gene stable ID',	'Protein stable ID',	'Gene name'])
rabbit_ortho_id_1 = rabbit_ortho_id_0.rename(columns = {'Gene stable ID': 'Rabbit_Gene_stable_ID', 'Protein stable ID':'Rabbit_Protein_stable_ID', 'Gene name':'Rabbit_Gene_name'})
rabbit_ortho_ids = rabbit_ortho_id_1.drop_duplicates() ## unique IDs of orthologous proteins

rabbit_ortho_ids.to_csv(Rabbit_file_with_ortho_ID, header=None)

#print(rabbit_ortho_ids)


with open(Rabbit_path_to_not_found_prot_ID, "w+") as not_found:
    wr = csv.writer(not_found)
    wr.writerow(rabbit_not_found)
#
############# Rabbit END ############
#
###### Dicty ######
#
dicty_not_found =[]
dicty_found_ortho_id = pd.DataFrame()
dicty_ortho_id_0 = pd.DataFrame()  ## For pheno_mapping

for i in Dicty_ortholog_prot:

    search_Dicty_db = df_Dicty.loc[(df_Dicty['Protein stable ID'] == i) ]

    if search_Dicty_db.empty:
        dicty_not_found.append(i)
#        print(" The following protein cannot be found in Database. Might be an old ID. sometimes there are -Retired ENSEMBL IDS-. search manually.Record entered in NOT FOUND file. \n\n"+i)
    else:
        dicty_found_ortho_id = dicty_found_ortho_id.append(search_Rabbit_db)

dicty_found_ortho_id.to_csv(path_to_Dicty_info_file)
dicty_ortho_id_0 = dicty_found_ortho_id.reindex(columns = ['Gene stable ID','Protein stable ID','Gene name'])
dicty_ortho_id_1 = dicty_ortho_id_0.rename(columns = {'Gene stable ID': 'Dicty_Gene_stable_ID', 'Protein stable ID':'Dicty_Protein_stable_ID', 'Gene name':'Dicty_Gene_name'})
dicty_ortho_ids = dicty_ortho_id_1.drop_duplicates() ## unique IDs of orthologous proteins

dicty_ortho_ids.to_csv(Dicty_file_with_ortho_ID, header=None)

#print(dicty_ortho_ids)

with open(Dicty_path_to_not_found_prot_ID, "w+") as not_found:
    wr = csv.writer(not_found)
    wr.writerow(dicty_not_found)
#
############# Dicty END ############
#

###### Rat ######
#
rat_not_found =[]
rat_found_ortho_id = pd.DataFrame()
rat_ortho_id_0 = pd.DataFrame()  ## For pheno_mapping

for i in Rat_ortholog_prot:

    search_Rat_db = df_Rat.loc[(df_Rat['Protein stable ID'] == i) ]

    if search_Rat_db.empty:
        rat_not_found.append(i)
    else:
        rat_found_ortho_id = rat_found_ortho_id.append(search_Rat_db)

rat_found_ortho_id.to_csv(path_to_Rat_info_file)
rat_ortho_id_0 = rat_found_ortho_id.reindex(columns = ['Gene stable ID','Protein stable ID','Gene name'])
rat_ortho_id_1 = rat_ortho_id_0.rename(columns = {'Gene stable ID': 'Rat_Gene_stable_ID', 'Protein stable ID':'Rat_Protein_stable_ID', 'Gene name':'Rat_Gene_name'})
rat_ortho_ids = rat_ortho_id_1.drop_duplicates() ## unique IDs of orthologous proteins

rat_ortho_ids.to_csv(Rat_file_with_ortho_ID, header=None)

#print(rat_ortho_ids)

with open(Rat_path_to_not_found_prot_ID, "w+") as not_found:
    wr = csv.writer(not_found)
    wr.writerow(rat_not_found)
#
############# Rat END ############
#
