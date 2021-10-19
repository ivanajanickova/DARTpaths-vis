#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 20:00:25 2020

@author: d
"""

##### DARTpaths 
##### script for data preparation for the orthology-mapping pipeline. 
#### First step after downloading Reactome file with Uniprot accesions.
#### This script is called from 00_prepare_data.sh 

### reference - websites
### https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/UniProt_programmatically_py3.pdf
### https://www.uniprot.org/help/api_idmapping
### https://www.uniprot.org/help/api
### https://www.uniprot.org/help/api_retrieve_entries
### https://www.uniprot.org/help/programmatic_access

### MAKE SURE YOU ARE CONNECTED TO INTERNET!

### This script takes 
## 1.file downloaded from Reactome
## 2. name of pathway

### using this, it will generate data that will be used by reconcile1.sh
### While it is running it retrieves data from Uniprot.
### makes necessary directories
### creates necessary files and keep in its location 

import os
import sys
import csv
import requests
from Bio import SeqIO

#######################  take this from bash/commandline

### specify name for pathway. use abbreviation if actual name is too long.
pathway_name = sys.argv[1]

## specified by 00_prepare_data.sh
path_to_pathway_files = sys.argv[2]
prots_file_name = sys.argv[3]

########################

file_complete = path_to_pathway_files+'/'+'Pathway_files/'+prots_file_name


list_of_prot_entry =[]

with open(file_complete,"r") as f:
    entry_line = f.read()
entry_0 = entry_line.strip()
entry = entry_0.split("\n")

#print(entry)
#print(len(entry))

## Exclude the header line
list_of_prot_entry = entry[1:]
#print(list_of_prot_entry)
#
#print(list_of_prot_entry[0].split('\t'))
#print(type(list_of_prot_entry[0]))
#list_of_accession =[]

#path_to_pathway_files = '/Users/d/Desktop/2020_0218_scripts/Mock_automatepipeline/DARTpaths/data_preparation/Pathways'
path_to_pathway = path_to_pathway_files + "/pathways/"+ pathway_name

###### prepare folders for keeping ortholog files.     #######################

Ortholog_folder = path_to_pathway + '/' + "Orthologs_" +pathway_name +"_pathway"

os.mkdir(Ortholog_folder)

path_to_pathway_ortholog_folder_root = path_to_pathway + '/' + "Orthologs_" +pathway_name +"_pathway/"

Mouse_folder = path_to_pathway_ortholog_folder_root + "Mouse_orthologs"
Rat_folder = path_to_pathway_ortholog_folder_root + "Rat_orthologs"
Zebrafish_folder = path_to_pathway_ortholog_folder_root +"Zebrafish_orthologs"
Rabbit_folder = path_to_pathway_ortholog_folder_root + "Rabbit_orthologs"
Cele_folder = path_to_pathway_ortholog_folder_root + "Cele_orthologs"
Dicty_folder = path_to_pathway_ortholog_folder_root + "Dicty_orthologs"


os.mkdir(Mouse_folder)
os.mkdir(Rat_folder)
os.mkdir(Zebrafish_folder)
os.mkdir(Rabbit_folder)
os.mkdir(Cele_folder)
os.mkdir(Dicty_folder)

#############################

fullURL=('http://www.uniprot.org/uniprot/?')

BASE = 'http://www.uniprot.org'
KB_ENDPOINT = '/uniprot/'
TOOL_ENDPOINT = '/uploadlists/'

def map_retrieve(ids2map, source_fmt='ACC+ID',target_fmt='STRING_ID', output_fmt='list'):
    if hasattr(ids2map, 'pop'):
        ids2map = ' '.join(ids2map)
    payload = {'from': source_fmt, 
               'to': target_fmt,
               'format': output_fmt,
               'query': ids2map,
               }
    response = requests.get(BASE + TOOL_ENDPOINT, params=payload)
    if response.ok:
        return response.text
    else:
        response.raise_for_status()
    
result = requests.get(fullURL)


for i in list_of_prot_entry:
    
    per_line = i
    single_entry = i.split('\t')
    molecule_type = single_entry[0]
    uniprot_acc = single_entry[1]
    uniprot_entry_gene = single_entry[2]
        
    string_id = map_retrieve(uniprot_acc, source_fmt='ACC')
    gene_name = map_retrieve(uniprot_acc,target_fmt='GENENAME')
    ENSEMBL_human_gene_id = map_retrieve(uniprot_acc,target_fmt='ENSEMBL_ID')
    
    pathway_prot_name = pathway_name+'_'+gene_name.strip()
    
    summary_of_this_prot=[]
    
    ### keep info of this protein in summary file
    summary_of_this_prot.append(pathway_prot_name)
    summary_of_this_prot.append(ENSEMBL_human_gene_id.strip())
    summary_of_this_prot.append(string_id.strip())
    
    ### get FASTA sequence
    uniprot_id = 'id:'+uniprot_acc

    payload = {'query': uniprot_id,'format': 'fasta'}
    uniprot_fastaseq = requests.get(BASE + KB_ENDPOINT, params=payload)
    
    
    path_to_pathway = path_to_pathway_files +'/'+'pathways'

    path_to_prot_folder = path_to_pathway+'/'+pathway_name+'/'+pathway_prot_name
    path_to_file_folder = path_to_prot_folder+'/'+pathway_prot_name+'_searches_literature'
    
    os.mkdir(path_to_prot_folder)
    os.mkdir(path_to_file_folder)
    
    ### kept under "PATHWAY_PROT_searches_literature" folder.
    ### files in which data from Uniprot will be stored. 
    prot_fasta_file = path_to_file_folder+'/'+pathway_prot_name+'_uniprot.fasta'
    prot_reactome_file = path_to_file_folder+'/'+pathway_prot_name+'_prots_reactome'
    prot_string_file = path_to_file_folder+'/'+pathway_prot_name+'_reactome.fa'
    
    prot_summary_file = path_to_file_folder+'/'+pathway_prot_name+'_summary.txt'
    
    with open(prot_summary_file, "w+") as summary_file:
        write_summary=csv.writer(summary_file)
        write_summary.writerow(summary_of_this_prot)
        
    with open(prot_reactome_file,"w+") as entry_file :
        entry_file.write(per_line[:])
    
    with open(prot_fasta_file,"w+") as f :
        f.write(uniprot_fastaseq.text[:])
    
    with open(prot_fasta_file,"r") as f, open(prot_string_file,"w+") as s :
        fasta_entry = SeqIO.parse(f,'fasta')
        for i in fasta_entry:
    #        print(i.description)
            i.id = string_id.strip()
            i.name = ""
            i.description = ""
        SeqIO.write(i,s,'fasta')   


    
