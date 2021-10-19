#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 17:07:22 2019

@author: d
"""

#### get_orthoinfo.py script re-written which uses dataframe from Pandas instead of dictionary.
####  much faster than dictionary, but might be even more efficient with sqlite OR TinyDB

import re
import sys
import time

import pandas as pd

protein_name=sys.argv[1]
current_path=sys.argv[2]
path_to_db_root = sys.argv[3]
query_protein_id = sys.argv[4]

prot = protein_name

line1 = "#{} protein family orthologs".format(prot)
line2 = "#{} Signature: Diksha Bhalla".format(time.strftime('%d-%m-%Y'))

#path_to_orthologfile = "/Users/d/Documents/PythonScripts/20191208/AHR_pathway_orthologs/"
full_path_to_ortholog_file = current_path+'/'+protein_name+"_orthologs"+'/'+protein_name+"_orthoinfo"+'/'+protein_name+"_human_orthologs.txt"
path_to_ortho_info_file = current_path+'/'+protein_name+"_orthologs"+'/'+protein_name+"_orthoinfo"+'/'+protein_name+"_QueryProt_orthoinfo.txt"

#query_protein_AHR1 = "9606.ENSP00000242057"
#query_protein_AHR2 = "9606.ENSP00000323816"
#query_protein_AHR3 = "9606.ENSP00000351407"
#query_protein_AHR4 = "9606.ENSP00000307479"
#query_protein_AHR5 = "9606.ENSP00000279146"
#query_protein_AHR6 = "9606.ENSP00000325875"
#query_protein_AHR7 = "9606.ENSP00000262033"

_ , query_protein_id_code = query_protein_id.split('.',1)

cele_path_to_db = path_to_db_root + "CaenorhabditisElegans_geneprot_database"
rabbit_path_to_db = path_to_db_root + "OryctolagusCuniculus_geneprot_database"
danio_path_to_db = path_to_db_root + "DanioRerio_geneprot_database"
rat_path_to_db = path_to_db_root + "RattusNorvegicus_geneprot_database"
mouse_path_to_db = path_to_db_root + "MusMusculus_geneprot_database"
dicty_path_to_db = path_to_db_root + "DictyosteliumDiscoideum_geneprot_database"
human_path_to_db = path_to_db_root + "HomoSapiens_geneprot_database"

df_cele = pd.read_csv(cele_path_to_db,delimiter='\t', header = 'infer')
df_Rabbit = pd.read_csv(rabbit_path_to_db,delimiter='\t', header = 'infer')
df_Danio = pd.read_csv(danio_path_to_db,delimiter='\t', header = 'infer')
df_Rat = pd.read_csv(rat_path_to_db,delimiter='\t', header = 'infer')
df_Mouse = pd.read_csv(mouse_path_to_db,delimiter='\t', header = 'infer')
df_Dicty = pd.read_csv(dicty_path_to_db,delimiter='\t', header = 'infer')
df_Human = pd.read_csv(human_path_to_db,delimiter='\t', header = 'infer')

with open(full_path_to_ortholog_file,'r') as ortho_tab:
	ortho_text = ortho_tab.read()

ortho_text.strip()
human_prots = ortho_text.split("\n\n")

list_of_query_protein = []
list_of_NON_query_protein = []

for i in human_prots:
    if i == "":
        continue
    else:
        human_prot_full, right_side = i.split('\t')
        human_prot_id, h_sp_name = human_prot_full.split('__')


        if human_prot_id == query_protein_id:
            list_of_query_protein.append(i)
        else:
            list_of_NON_query_protein.append(i)

search_Mus_db = df_Mouse.loc[(df_Mouse['homology_species'] == "homo_sapiens") & (df_Mouse['homology_protein_stable_id'] == query_protein_id_code)]
search_cele_db = df_cele.loc[(df_cele['homology_species'] == "homo_sapiens") & (df_cele['homology_protein_stable_id'] == query_protein_id_code)]
search_Dicty_db = df_Dicty.loc[(df_Dicty['homology_species'] == "homo_sapiens") & (df_Dicty['homology_protein_stable_id'] == query_protein_id_code)]
search_Danio_db = df_Danio.loc[(df_Danio['homology_species'] == "homo_sapiens") & (df_Danio['homology_protein_stable_id'] == query_protein_id_code)]
search_Human_db = df_Human.loc[(df_Human['homology_species'] == "homo_sapiens") & (df_Human['homology_protein_stable_id'] == query_protein_id_code)]


search_Rabbit_db = df_Rabbit.loc[(df_Rabbit['Protein stable ID'] == query_protein_id_code)]
search_Rat_db = df_Rat.loc[(df_Rat['Protein stable ID'] == query_protein_id_code)]

frames_001 =[search_Human_db,search_Mus_db,search_cele_db, search_Dicty_db,search_Danio_db ]
merged_frame_001 = pd.concat(frames_001)

## combine Rat and Rabbit  ## same format
frames_002 = [search_Rabbit_db,search_Rat_db]

with open(path_to_ortho_info_file,"w+") as f:
    f.write(line1 + "\n" + line2 + '\n\n')

merged_frame_001.to_csv(path_to_ortho_info_file,mode = 'a',columns=['homology_gene_stable_id','homology_protein_stable_id','homology_species',  'homology_type' ,'species','gene_stable_id','protein_stable_id','homology_identity'])


with open(path_to_ortho_info_file,"a") as f_1:
    f_1.write('\n\n')
    search_Rabbit_db.to_csv(path_to_ortho_info_file,mode = 'a',columns=['Gene stable ID','Protein stable ID','Rabbit homology type','Rabbit gene stable ID','Rabbit protein or transcript stable ID','Rabbit gene name','%id. query gene identical to target Rabbit gene','Rabbit orthology confidence [0 low, 1 high]'])
    f_1.write('\n\n')
    search_Rat_db.to_csv(path_to_ortho_info_file,mode = 'a',columns=['Gene stable ID','Protein stable ID','Rat homology type','Rat gene stable ID','Rat protein or transcript stable ID','Rat gene name','%id. query gene identical to target Rat gene','Rat orthology confidence [0 low, 1 high]'])
